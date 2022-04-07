function gather_useful_info

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
all_soz_locs = cell(npts,1);
all_soz_lats = cell(npts,1);
all_soz_bin = cell(npts,1);
all_spikes = cell(npts,1);
all_rl = cell(npts,1);
all_ns = cell(npts,1);
all_labels = cell(npts,1);
all_names = cell(npts,1);
all_fc = cell(npts,1);
all_locs = cell(npts,1);
all_coh = cell(npts,1);

%% Loop over patients
for p = 1:npts

    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    ns =summ.ns_car;
    spikes = summ.spikes;
    rl = summ.rl;
    labels = summ.labels;
    soz_loc = summ.soz.loc;
    soz_lat = summ.soz.lat;
    soz_labels = summ.soz.labels;
    name = summ.name;
    fc = summ.avg_fc;
    locs = summ.locs;
    coh = wrap_or__unwrap_adjacency_fc_toolbox(summ.avg_coh);
    
    all_names{p} = name;
    
    %% Get soz
    all_soz_locs{p} = soz_loc;
    all_soz_lats{p} = soz_lat;
    
    % binary SOZ identification
    soz_bin = zeros(length(labels),1);
    soz_bin(ismember(labels,soz_labels)) = 1;
    
    %% Remove non-intracranial
    ekg = find_non_intracranial(labels);
    ns(ekg,:) = [];
    spikes(ekg,:) = [];
    labels(ekg) = [];
    soz_bin(ekg) = [];
    rl(ekg,:) = [];
    fc(ekg,:) = []; fc(:,ekg) = [];
    locs(ekg,:) = [];
    coh(ekg,:,:) = [];
    coh(:,ekg,:) = [];
    
    %% SOZ bin
    all_soz_bin{p} = soz_bin;
    all_labels{p} = labels;
    
    %% Get averages over times
    avg_spikes = nanmean(spikes,2);
    avg_ns = nanmean(ns,2);
    avg_rl = nanmean(rl,2);
    
    all_spikes{p} = avg_spikes;
    all_ns{p} = avg_ns;
    all_rl{p} = avg_rl;
    all_fc{p} = fc;
    all_locs{p} = locs;
    all_coh{p} = coh;
    
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
out.all_soz_locs = all_soz_locs;
out.all_soz_lats = all_soz_lats;
out.all_soz_bin = all_soz_bin;
out.all_spikes = all_spikes;
out.all_locs = all_locs;
out.all_rl = all_rl;
out.all_fc = all_fc;
out.all_ns = all_ns;
out.all_coh = all_coh;
out.all_labels = all_labels;
out.all_names = all_names;


%% Save
save([out_folder,'main_out.mat'],'out');
    
end