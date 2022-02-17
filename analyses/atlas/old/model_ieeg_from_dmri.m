function model_ieeg_from_dmri

%{
Big idea: 
- ieeg atlas problematic due to sparsity of data
- fmri and structural connectivity reasonably correlated with ieeg
connectivity
- can i build a model to predict ieeg edge-level connectivity based on this (less
sparse data)?
- and then compare the deviation between actual patient edge-level
connectivity and the predicted connectivity and see if this deviation
localizes the SOZ?
%}


which_atlas = 'brainnetome';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];
bct_folder= locations.bct;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load ieeg atlas
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;

atlas = out.atlas;
names = out.atlas_names;
pt_names = out.pt_names;

if size(names,1) < size(names,2)
    names = names';
end
ieeg = out.atlas;
all_ieeg = ieeg;
ieeg = nanmean(ieeg,3);
sozs = out.sozs;

%% convert sozs to binary matrix
nregions = size(ieeg,1);
all_regions = 1:nregions;
bin_soz = cell2mat(cellfun(@(x) ismember(all_regions,x),sozs,'uniformoutput',false));

if 0
    imagesc(bin_soz)
    xticks(1:size(bin_soz,2))
    yticks(1:size(bin_soz,1))
    yticklabels(pt_names)
    xticklabels(names)
end

%% Load fmri atlases
fmri = load([atlas_folder,which_atlas,'_sc_fc.mat']);
fmri = fmri.out;
fmri_fc = (fmri.FC);
fmri_sc = fmri.SC;

%% Find and remove nan rows and things not usually targeted with implantation
basal_ganglia = contains(names,'BG_');
thalamus = contains(names,'Tha_');
nan_rows = (sum(isnan(ieeg),1) == size(ieeg,2))';
remove = basal_ganglia | thalamus | nan_rows;

fmri_fc = fmri_fc(~remove,~remove);
fmri_sc = fmri_sc(~remove,~remove);
ieeg = ieeg(~remove,~remove);
all_ieeg = all_ieeg(~remove,~remove,:);
bin_soz =bin_soz(:,~remove);
names = names(~remove);

%% Flatten
ieegw = wrap_or_unwrap_adjacency_fc_toolbox(ieeg);
fmri_fcw = wrap_or_unwrap_adjacency_fc_toolbox(fmri_fc);
fmri_scw = wrap_or_unwrap_adjacency_fc_toolbox(fmri_sc);

%% Measure sparsity
all_ieegw = wrap_or_unwrap_adjacency_fc_toolbox(all_ieeg);
missing_rows = squeeze(sum(isnan(all_ieeg),1) == size(all_ieeg,1));
npts = size(all_ieeg,3);
nedges = size(all_ieegw,1);

pts_per_edge = measure_sparsity(all_ieegw);
if 0
    histogram(pts_per_edge)
end



%% Linear model
T = table(ieegw,fmri_scw,fmri_fcw);
lm = fitlm(T,'ieegw ~ fmri_scw + fmri_fcw');
params = lm.CoefficientNames;

%% Generate prediction

prediction = lm.Coefficients.Estimate(1)*ones(nedges,1) + ...
    lm.Coefficients.Estimate(2)*fmri_scw + ...
    lm.Coefficients.Estimate(3)*fmri_fcw;
prediction_uw = wrap_or_unwrap_adjacency_fc_toolbox(prediction);

if 0
    figure
    imagesc(prediction_uw)
    r = corr(prediction,ieegw,'rows','pairwise');
end

%% Measure residuals for each patient
square_resid = nan(npts,nedges);
for ip = 1:npts
    curr_ieeg = all_ieegw(:,ip);
    resid = curr_ieeg-prediction;
    square_resid(ip,:) = (resid);
end

%% Unwrap
square_resid_uw = wrap_or_unwrap_adjacency_fc_toolbox(square_resid');

%% Add nans back in where appropriate
for ip = 1:npts
    square_resid_uw(missing_rows(:,ip),:,ip) = nan;
    square_resid_uw(:,missing_rows(:,ip),ip) = nan;
end

%% Get mean residual across connections
ns_resid = squeeze(nanmean(square_resid_uw,1));
ns_resid = ns_resid';

%% Some visualization


if 0
    figure
    skip = 5;
    turn_nans_gray(ns_resid)
    yticks(1:size(ns_resid,1))
    xticks(1:skip:size(ns_resid,2))
    xticklabels(names(1:skip:length(names)))
    yticklabels(pt_names)
end

if 1
    figure
    for i = 1:npts
        curr_resid = ns_resid(i,:);
        curr_soz = bin_soz(i,:);
        nan_elements = isnan(curr_resid);
        curr_resid(nan_elements) = [];
        curr_soz(nan_elements) = [];
        [curr_resid,I] = sort(curr_resid);
        curr_soz = curr_soz(I);
        plot(1,1)
        hold off
        plot(curr_resid,'ko');
        hold on
        plot(find(curr_soz),curr_resid(curr_soz),'ro')
        xticks(1:length(curr_soz))
        xticklabels(names(~nan_elements))
        %names(~nan_elements)
        
        pause
        hold off
    end
end


%% Plot orders
plot_orders_mats(ns_resid',bin_soz')


if 0
    figure
    imagesc(square_resid)
end



end