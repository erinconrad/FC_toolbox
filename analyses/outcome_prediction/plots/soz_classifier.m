%{

To do:
- add ability to use only sleep

%}

function [T_test,T_train,all_auc] = soz_classifier
nb = 10;
params.only_sleep = 0;
params.models = {'ana','ana_cov','ana_cov_spikes','ana_cov_spikes_ns','ana_cov_ns'};
params.pretty_name = {'Anatomy','Add coverage density','Add spike rates','Add connectivity','Add connectivity'};
params.do_norm = 1;
params.which_atlas = 'aal_bernabei';
params.max_spikes = 0.1;
params.prop_train = 2/3;

%% File locations
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];

nmodels = length(params.models);
all_auc = nan(nb,nmodels);

for im = 1:nmodels
    model_info(im).name = params.models{im};
    model_info(im).pretty_name = params.pretty_name{im};
    model_info(im).all_auc = nan(nb,1);
    model_info(im).all_roc = cell(nb,1);
end

for ib = 1:nb
    if mod(ib,10) == 0, fprintf('\nDoing %d of %d\n',ib,nb); end
    out = individual_classifier(params);
    for im = 1:nmodels
        all_auc(ib,im) = out.model(im).auc;
        model_info(im).all_auc(ib) = out.model(im).auc;
        model_info(im).all_roc{ib} = out.model(im).roc;
    end
end

%array2table(all_auc,'variablenames',models)
for im = 1:nmodels
    [unify_x,unify_y] = unify_roc(model_info(im).all_roc);
    model_info(im).x = unify_x;
    model_info(im).y = unify_y;
    model_info(im).ym = median(unify_y);
    model_info(im).ly = prctile(unify_y,25,1);
    model_info(im).uy = prctile(unify_y,75,1);
    
end

figure
set(gcf,'position',[10 10 1400 350])
tiledlayout(1,5,'tilespacing','tight','padding','tight')

for im = [1 2 5 3 4]
    nexttile
    mp = shaded_error_bars_fc(model_info(im).x,model_info(im).ym,...
        [model_info(im).ly',model_info(im).uy'],'k');
    hold on
    plot([0 1],[0 1],'k--','linewidth',2)
    title(model_info(im).pretty_name)
    legend(mp,sprintf('Median AUC = %1.2f',median( model_info(im).all_auc)),...
        'location','southeast','fontsize',15)
    set(gca,'fontsize',15)
    xlabel('False positive rate')
    ylabel('True positive rate')
end

print(gcf,[plot_folder,'prediction_nospikes'],'-dpng')

end

function out = individual_classifier(params)
%{
Predict whether an electrode is in the SOZ based on 1) its spike rate
normalized across other electrodes within the patient and 2) its anatomical
localization (a categorical variable that can be MT, TN, WM, or OC)
%}



%% Parameters
max_spikes = params.max_spikes; % max spikes/elecs/min to include in model
which_atlas = params.which_atlas;
prop_train = params.prop_train;
do_norm = params.do_norm;
models = params.models;


%% File locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
spike_folder = [results_folder,'analysis/outcome/data/'];
atlas_folder = [results_folder,'analysis/atlas/'];


%% Load out file and get roc stuff
out = load([spike_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
soz = out.all_soz_bin;
npts = length(soz);
labels = out.all_labels;
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
atlas_nums = atlas_out.atlas_nums;
atlas_names = atlas_out.atlas_names;

%% Spatially normalize the FC matrix
%[resid,f] = fit_distance_model(locs,fc,soz,rate,max_spikes);
%ns_resid = cellfun(@(x) nanmean(x,2),resid,'uniformoutput',false);
%ns = ns_resid;
ns = cellfun(@(x) nanmean(x,2),fc,'uniformoutput',false);

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



%% Get stuff into friendly format
vec_rate = [];
vec_loc = {};
vec_soz = [];
vec_pt_idx = [];
vec_ns = [];
%vec_pred_conn_D = [];
vec_dens = [];

for ip = 1:npts
    curr_rate = rate{ip};
    curr_soz = soz{ip};
    curr_ana_loc = elec_broad{ip};
    curr_ns = ns{ip};
    curr_loc = locs{ip};
    
    % predict connectivty just based on locs (using model from above)
    %{
    predict_conn = predict_conn_from_loc(f,curr_loc);
    predict_conn = nanmean(predict_conn,2);
    vec_pred_conn_D = [vec_pred_conn_D;predict_conn];
    %}
    
    % get density
    density = estimate_coverage_density(curr_loc,1e8);
        
    % normalize within patient
    if do_norm
        curr_rate = (curr_rate-nanmean(curr_rate))./nanstd(curr_rate);
        curr_ns = (curr_ns-nanmean(curr_ns))./nanstd(curr_ns);
        density = (density - nanmean(density))./nanstd(density);
    end
    
    vec_rate = [vec_rate;curr_rate];
    vec_dens = [vec_dens;density];
    vec_loc = [vec_loc;curr_ana_loc];
    vec_soz = [vec_soz;curr_soz];
    vec_ns = [vec_ns;curr_ns];
    vec_pt_idx = [vec_pt_idx;repmat(ip,length(curr_soz),1)];
end

%% Make table
T = table(vec_soz,vec_pt_idx,vec_rate,vec_loc,vec_ns,vec_dens);

%% Remove nan and inf rows
nan_rows = isnan(T.vec_soz) | isnan(T.vec_pt_idx) | isnan(T.vec_rate) | ...
    isnan(T.vec_ns) | isnan(T.vec_dens);
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

for im = 1:length(models)
    curr_model = models{im};
    switch curr_model
        
        case 'ana'
            formula = 'vec_soz ~ vec_loc';
            formula_me = 'vec_soz ~ vec_loc + vec_pt_idx';
        case 'ana_cov'
            formula = 'vec_soz ~ vec_loc + vec_dens';
            formula_me = 'vec_soz ~ vec_loc + vec_dens+ vec_pt_idx';
        case 'ana_cov_spikes'
            formula = 'vec_soz ~ vec_loc + vec_dens + vec_rate';
            formula_me = 'vec_soz ~ vec_loc + vec_dens + vec_rate + vec_pt_idx';
        case 'ana_cov_spikes_ns'
            formula = 'vec_soz ~ vec_loc + vec_dens + vec_rate + vec_ns';
            formula_me = 'vec_soz ~ vec_loc + vec_dens + vec_rate + vec_ns + vec_pt_idx';
        case 'ana_cov_ns'
            formula = 'vec_soz ~ vec_loc + vec_dens + vec_ns';
            formula_me = 'vec_soz ~ vec_loc + vec_dens + vec_ns + vec_pt_idx';
    end
    
    %% Do both the glm and glme
    glm = fitglm(T_train,formula,'Distribution','Poisson','Link','log');
    glme = fitglme(T_train,formula_me,'Distribution','Poisson','Link','log');
    
    %% Test glm on testing data
    params = glm.CoefficientNames;
    
    % Calculate the model response
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
    
    % Compare classification against true
    all_soz = T_test.vec_soz==1;
    all_no_soz = T_test.vec_soz==0;
    [X,Y,T,AUC,opt] = perfcurve(T_test.vec_soz,classification,1);
    
    out.model(im).name = curr_model;
    out.model(im).glme = glme;
    out.model(im).roc = [X,Y];
    out.model(im).auc = AUC;
    
end





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
