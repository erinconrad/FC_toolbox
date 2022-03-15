function corr_spikes_fc

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
int_folder = [results_folder,'analysis/intermediate/'];
out_folder = [results_folder,'analysis/outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% Initialize result variables
all_corrs_sp = cell(npts,1);
all_corrs_pear = cell(npts,1);
all_avg_corrs_sp = nan(npts,1);
all_avg_corrs_pear = nan(npts,1);

%% Loop over patients
for p = 1:npts

    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    ns =summ.ns_car;
    spikes = summ.spikes;
    labels = summ.labels;
    
    %% Remove non-intracranial
    ekg = find_non_intracranial(labels);
    ns(ekg,:) = [];
    spikes(ekg,:) = [];
    labels(ekg) = [];
    
    %% Get averages over times
    avg_spikes = nanmean(spikes,2);
    avg_ns = nanmean(ns,2);
    
    %% Correlate average ns and spikes
    avg_corr_sp = corr(avg_spikes,avg_ns,'rows','pairwise','type','spearman');
    avg_corr_pear = corr(avg_spikes,avg_ns,'rows','pairwise','type','pearson');
    all_avg_corrs_sp(p) = avg_corr_sp;
    all_avg_corrs_pear(p) = avg_corr_pear;
    
    %% Correlate time varying
    ntimes = size(spikes,2);
    corr_sp = nan(ntimes,1);
    corr_pear = nan(ntimes,1);
    for it = 1:ntimes
        corr_sp(it) = corr(spikes(:,it),ns(:,it),'rows','pairwise','type','spearman');
        corr_pear(it) = corr(spikes(:,it),ns(:,it),'rows','pairwise','type','pearson');
    end
    
    all_corrs_sp{p} = corr_sp;
    all_corrs_pear{p} = corr_pear;

end

out.avg_corr_sp = all_avg_corrs_sp;
out.avg_corr_pear = all_avg_corrs_pear;
out.all_corrs_sp = all_corrs_sp;
out.all_corrs_pear = all_corrs_pear;

%% Save
save([out_folder,'corr_out.mat'],'out');
    
end