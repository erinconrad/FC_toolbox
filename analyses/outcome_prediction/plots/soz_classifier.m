function [T_test,T_train,all_auc] = soz_classifier
nb = 1;
all_auc = nan(nb,1);


for ib = 1:nb
    soz_roc_out = individual_classifier;
    all_auc(ib) = soz_roc_out.auc;
end

T_test = soz_roc_out.T_test;
T_train = soz_roc_out.T_train;
all_auc

end

function soz_roc_out = individual_classifier
%{
Predict whether an electrode is in the SOZ based on 1) its spike rate
normalized across other electrodes within the patient and 2) its anatomical
localization (a categorical variable that can be MT, TN, WM, or OC)
%}

%% Parameters
which_atlas = 'brainnetome';%'aal_bernabei';
broad_locs = {'mesial temporal','temporal neocortical','other cortex'};
prop_train = 2/3;
do_norm = 1; % normalize within patient
do_anat = 0;


%% File locations
locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
spike_folder = [results_folder,'analysis/outcome/data/'];
atlas_folder = [results_folder,'analysis/atlas/'];


%% Load out file and get roc stuff
out = load([spike_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
rl = out.all_rl;
soz = out.all_soz_bin;
loc = out.all_soz_locs;
npts = length(soz);
labels = out.all_labels;
ns = out.all_ns;
fc = out.all_fc;
locs = out.all_locs;

%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);



%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

%% get atlas stuff
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;
spikes_atlas = atlas_out.spikes_atlas;
atlas_nums = atlas_out.atlas_nums;
atlas_names = atlas_out.atlas_names;

%% Spatially normalize the FC matrix
resid = fit_distance_model(locs,fc,soz);
ns_resid = cellfun(@(x) nanmean(x,2),resid,'uniformoutput',false);

%% Localize regions into broad categories
broad = localize_regions(atlas_names,which_atlas);
broad_no_lat = cell(length(broad),1);
for ib = 1:length(broad)
    curr = broad{ib};
    if isempty(curr),continue; end
    curr = strrep(curr,'left ','');
    curr = strrep(curr,'right ','');
    broad_no_lat{ib} = curr;
end

%% Get localizations for each patient
elec_broad = cell(npts,1);
for ip = 1:npts
    elec_broad{ip} = get_anatomy_elecs(labels{ip},atlas_elec_labels{ip},atlas_elec_regions{ip},broad_no_lat,atlas_nums);
end

%% Normalize spike rates by anatomical location
z = cell(npts,1);
for ip = 1:npts
    z{ip} = normalize_spike_rates(labels{ip},atlas_elec_labels{ip},...
        atlas_elec_regions{ip},spikes_atlas,rate{ip},atlas_nums);
end
rate_z = z;


%% Get stuff into friendly format
vec_rate = [];
vec_loc = {};
vec_soz = [];
vec_pt_idx = [];
vec_ns = [];
vec_rl = [];
vec_ns_resid = [];
vec_rate_z = [];

for ip = 1:npts
    curr_rate = rate{ip};
    curr_soz = soz{ip};
    curr_loc = elec_broad{ip};
    curr_ns = ns{ip};
    curr_rl = rl{ip};
    curr_ns_resid = ns_resid{ip};
    curr_rate_z = rate_z{ip};
        
    % normalize within patient
    if do_norm
        curr_rate = (curr_rate-nanmean(curr_rate))./nanstd(curr_rate);
        curr_ns_resid = (curr_ns_resid - nanmean(curr_ns_resid))./nanstd(curr_ns_resid);
        curr_rate_z = (curr_rate_z-nanmean(curr_rate_z))./nanstd(curr_rate_z);
        curr_rl = (curr_rl-nanmean(curr_rl))./nanstd(curr_rl);
    end
    
    vec_rate = [vec_rate;curr_rate];
    vec_loc = [vec_loc;curr_loc];
    vec_soz = [vec_soz;curr_soz];
    vec_rl = [vec_rl;curr_rl];
    vec_ns = [vec_ns;curr_ns];
    vec_rate_z = [vec_rate_z;curr_rate_z];
    vec_ns_resid = [vec_ns_resid;curr_ns_resid];
    vec_pt_idx = [vec_pt_idx;repmat(ip,length(curr_soz),1)];
end

%% Make table
T = table(vec_soz,vec_pt_idx,vec_rate,vec_loc,vec_ns,vec_rl,vec_ns_resid,vec_rate_z);

%% Remove nan and inf rows
nan_rows = isnan(T.vec_soz) | isnan(T.vec_pt_idx) | isnan(T.vec_rate) | ...
    isnan(T.vec_rl) | isnan(T.vec_ns) | isnan(T.vec_ns_resid);
T(nan_rows,:) = [];

%% Remove rows without localization
no_loc = cellfun(@isempty,T.vec_loc);
T(no_loc,:) = [];

%% Divide into training and testing set
training = randsample(npts,floor(npts*prop_train));
training_idx = ismember(T.vec_pt_idx,training);
testing_idx = ~ismember(T.vec_pt_idx,training);

T_train = T(training_idx,:);
T_test = T(testing_idx,:);

%% Confirm that I separated by patients
train_pts = unique(T_train.vec_pt_idx);
test_pts = unique(T_test.vec_pt_idx);
assert(isempty(intersect(train_pts,test_pts)))

if do_anat
glm = fitglm(T_train,...
    'vec_soz ~ vec_rate + vec_ns',...
    'Distribution','Poisson','Link','log');
else

glm = fitglm(T_train,...
    'vec_soz ~ vec_rate',...
    'Distribution','Poisson','Link','log');
end

%% Test model on testing data
% Find the reference category
params = glm.CoefficientNames;

%{
non_ref_cats = {};
for p  = 2:length(params)
    if ~contains((params{p}),'vec_loc')
        continue
    end
    
    C = strsplit((params{p}),'_');
    non_ref_cats = [non_ref_cats;C{3}];
    
end
%}


asum = zeros(size(T_test,1),1);
asum = asum + glm.Coefficients.Estimate(1);
for p = 2:length(params)
    est = glm.Coefficients.Estimate(p);
    if contains((params{p}),'vec_loc')
        C = strsplit((params{p}),'_');
        cat = C{3};
        
        % find the indices of T_test that fit this category
        match = strcmp(T_test.vec_loc,cat);
        asum(match) = asum(match) + est;
    else
        
        asum = asum +  T_test.(params{p})*est;
    end
end


classification = logistic(asum);
all_soz = T_test.vec_soz==1;
all_no_soz = T_test.vec_soz==0;
class_soz = classification(all_soz);
class_no_soz = classification(all_no_soz);
[roc,auc,disc,disc_I] = calculate_roc(class_no_soz,class_soz,1e3);

soz_roc_out.roc = roc;
soz_roc_out.auc = auc;
soz_roc_out.T_train = T_train;
soz_roc_out.T_test = T_test;
soz_roc_out.glm = glm;

end

function elec_broad = get_anatomy_elecs(sp_labels,atlas_labels,regions,broad_no_lat,atlas_nums)
    
%% Remove ekg from atlas
ekg = find_non_intracranial(atlas_labels);
atlas_labels(ekg) = [];
regions(ekg) = [];

%% Reconcile labels
assert(isequal(sp_labels,atlas_labels)) 

elec_broad = cell(length(sp_labels),1);

% Loop over electrodes
for ich = 1:length(sp_labels)
    
    % get the region of this electrode
    curr_region = regions(ich);
    
    % get the index of this
    curr_idx = find(atlas_nums == curr_region);
    
    if isempty(curr_idx), continue; end
    
    elec_broad{ich} = broad_no_lat{curr_idx};
end

if 0
    table(sp_labels,elec_broad)
end


end

function out = logistic(x)

out = 1./(1+exp(-x));

end

function z = normalize_spike_rates(sp_labels,atlas_labels,regions,atlas,spikes,atlas_nums)
    
%% Remove ekg from atlas
ekg = find_non_intracranial(atlas_labels);
atlas_labels(ekg) = [];
regions(ekg) = [];

%% Reconcile labels
assert(isequal(sp_labels,atlas_labels)) 

z = nan(length(spikes),1);

% Loop over electrodes
for ich = 1:length(spikes)
    
    % get the region of this electrode
    curr_region = regions(ich);
    
    % get the index of this
    curr_idx = find(atlas_nums == curr_region);
    
    if isempty(curr_idx), z(ich) = nan; continue; end
    
    % normalize by atlas (across all patients for that region)
    z(ich) = (spikes(ich)-nanmedian(atlas(curr_idx,:)))./iqr(atlas(curr_idx,:));
    %z(ich) = (spikes(ich)-nanmean(atlas(curr_idx,:)))./nanstd(atlas(curr_idx,:));
end

if 0
    table(sp_labels,spikes,z)
end

end
