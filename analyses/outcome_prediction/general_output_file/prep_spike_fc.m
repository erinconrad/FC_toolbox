function prep_spike_fc

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [results_folder,'analysis/outcome/data/spike_corrs/'];
out_folder = [results_folder,'analysis/outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load main out file
main_out = load([out_folder,'main_out.mat']);
main_out = main_out.out;
npts = length(main_out.all_names);

%% Initialize spike info
sp_names = main_out.all_names;
all_chs_corr = cell(npts,1);
sp_chs_corr = cell(npts,1);

for ip = 1:npts
    
    name = sp_names{ip};
    
    %% Load the file
    if exist([data_folder,name,'.mat'],'file') == 0
        continue;
    end
    spikes = load([data_folder,name,'.mat']);
    spikes = spikes.spout;
    
    all_chs_corr{ip} = spikes.all_chs_corr;
    sp_chs_corr{ip} = spikes.sp_chs_corr;
    
end

%% Save file
spikes_out.names = sp_names;
spikes_out.all_chs_corr = all_chs_corr;
spikes_out.sp_chs_corr = sp_chs_corr;
save([out_folder,'spikes_out.mat'],'spikes_out');

end