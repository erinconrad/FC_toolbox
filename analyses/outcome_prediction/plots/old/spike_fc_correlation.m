function spike_fc_correlation

%{
Calculate the correlation between spikes and FC and the pre-spike vs
during-spike FC
%}

%% Locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))

%% Load out file
out = load([locations.paper_data_folder,'main_out.mat']);
out = out.out;

%% Also load spike out file
spout = load([locations.paper_data_folder,'spikes_out.mat']);
spout = spout.spikes_out;

%% Get total n for main result n
npts = length(out.all_locs);
any_locs = zeros(npts,1);
for ip = 1:npts
   if any(~isnan(out.all_locs{ip}),'all')
       any_locs(ip) = 1;
   end
end
nout.pts_with_any_locs = (any_locs);
nout.names = out.all_names;

all_chs_corr = spout.all_chs_corr;
sp_chs_corr = spout.sp_chs_corr;
single_ch_corr = spout.single_ch_corr;

%% Get stuff
rate = out.all_spikes;
good_spikes = out.good_spikes;
soz = out.all_soz_bin;
npts = length(soz);
labels = out.all_labels;
fc = out.all_fc;

%% Do individual spike analysis
all_sp_chs_corr = nan(npts,2); % first column is pre spike and 2nd column is spike
all_all_chs_corr = nan(npts,2);
all_single_chs_corr = nan(npts,2);
for ip = 1:npts
    curr_all = all_chs_corr{ip};
    curr_sp = sp_chs_corr{ip};
    curr_single = single_ch_corr{ip};
    
    if isempty(curr_all), continue; end
    
    % take mean across spikes
    curr_all_mean = nanmean(curr_all,1);
    curr_sp_mean = nanmean(curr_sp,1);
    curr_single_mean = nanmean(curr_single,1);
    
    all_sp_chs_corr(ip,:) =curr_sp_mean;
    all_all_chs_corr(ip,:) = curr_all_mean;
    all_single_chs_corr(ip,:) = curr_single_mean;
end
nout.all_sp_chs_corr = all_sp_chs_corr;
nout.all_all_chs_corr = all_all_chs_corr;
nout.all_single_chs_corr = all_single_chs_corr;

%% Get spike-fc correlations
spike_fc_corr = nan(npts,1);
for ip = 1:npts
    curr_labels = labels{ip};
    assert(sum(find_non_intracranial(curr_labels)) == 0)
    
    curr_avg_fc = nanmean(fc{ip},2);
    curr_spikes = rate{ip};
    spike_fc_corr(ip) = corr(curr_avg_fc,curr_spikes,'rows','pairwise','type','spearman');
    
end

nout.spike_fc_corr = spike_fc_corr;
nout.good_spikes = good_spikes;

save([locations.paper_plot_folder,'spike_analysis.mat'],'nout');

end