function investigate_re

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
subplot_path = [plot_folder,'ai_subplots/'];
if ~exist(subplot_path,'dir')
    mkdir(subplot_path)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
mt_data = load([inter_folder,'mt_out.mat']);
mt_data = mt_data.out;

which_sleep_stage = 1; % all = 1, wake =2, sleep = 3;
which_montage = 2;
   

spikes = mt_data.all_spikes(:,which_montage,which_sleep_stage);
re = mt_data.all_re(:,which_montage,which_sleep_stage);
npts = length(spikes);
spikes_re_corr = nan(npts,1);
for i = 1:npts
    if isempty(re{i}), continue; end
    re_4 = re{i}(:,:,4);
    re_4_ns = nanmean(re_4,2);
    spikes_curr = spikes{i};
    spikes_re_corr(i) = corr(spikes_curr,re_4_ns,'rows','pairwise');
    
end

        

end