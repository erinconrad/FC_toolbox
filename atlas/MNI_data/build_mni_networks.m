function build_mni_networks

%% Parameters
tw = 2; % 2 second time window for pc calculations
broad_regions = {'left mesial temporal','right mesial temporal',...
    'left temporal neocortical','right temporal neocortical',...
    'left other cortex','right other cortex'};
nb = length(broad_regions);

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/atlas/mni/'];
int_folder = [results_folder,'analysis/intermediate/'];
data_folder = [locations.main_folder,'data/'];
mni_folder = [data_folder,'MNI_open_ieeg/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load mni data
mni = load([mni_folder,'MatlabFile.mat']);

%% Get the eeg data
eeg = mni.Data_W;
fs = mni.SamplingFrequency;
hemisphere = mni.Hemisphere;
channelRegion = mni.ChannelRegion;
region_name = mni.RegionName;

% convert hemisphere to left/right
new_hemisphere = hemisphere;
for i = 1:length(hemisphere)
    if strcmp(hemisphere{i},'L')
        new_hemisphere{i} = 'left';
    elseif strcmp(hemisphere{i},'R')
        new_hemisphere{i} = 'right';
    end
end
hemisphere = new_hemisphere;

%% find conversion between MNI regions and my broad regions
regions = mniregions_to_region(region_name);

%% Loop over patients
% find unique patients
pts = unique(mni.Patient);
npts = length(pts);
pt_nets = cell(npts,1);
pt_regions = cell(npts,1);
pt_labels = cell(npts,1);
broad_connectivity = nan(nb,npts);
for i = 1:npts
    ip = pts(i); % select current patient
    
    % find corresponding channels
    curr_idx = mni.Patient == ip;
    
    % get current eeg, labels, hemispheres, regions
    curr_eeg = eeg(:,curr_idx);
    old_eeg = curr_eeg;
    curr_labels = mni.ChannelName(curr_idx);
    curr_region_num = channelRegion(curr_idx);
    curr_hemisphere = hemisphere(curr_idx);
    
    if length(curr_labels) == 2, continue; end
    
    pt_labels{i} = curr_labels;
    
    % do a CAR
    curr_eeg = curr_eeg - nanmean(curr_eeg,2);
    
    if 0
        show_eeg(curr_eeg,fs,curr_labels)
    end
    
    % Get PC network
    curr_net = pc_vector_calc(curr_eeg,fs,tw);
    
    % unwrap
    curr_net = wrap_or_unwrap_adjacency_fc_toolbox(curr_net);
    
    % fill up cell
    pt_nets{i} = curr_net;
    
    % Get full broad region/hemispheric combo of each channel
    curr_region = regions(curr_region_num);
    curr_region = cellfun(@(x,y) [x,' ',y],curr_hemisphere,curr_region,'Uniformoutput',false);
    pt_regions{i} = curr_region;
    
    % get the intraregional connectivity
    for ib = 1:nb
        curr_in_region = strcmp(curr_region,broad_regions{ib});
        broad_connectivity(ib,i) = nanmean(curr_net(curr_in_region,curr_in_region),'all');
    end
    
end

out.pt_ids =pts;
out.pt_labels = pt_labels;
out.pt_regions = pt_regions;
out.pt_nets = pt_nets;
out.broad_connectivity = broad_connectivity;

if 0
    figure
    turn_nans_gray(broad_connectivity)
    yticks(1:6)
    yticklabels(broad_regions)
end

save([out_folder,'mni_net.mat'],'out');

end