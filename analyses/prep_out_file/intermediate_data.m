function intermediate_data

%% Parameters
m = 2; % do not change

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/intermediate/'];
spikes_folder = [results_folder,'all_out/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

validation_file = [scripts_folder,'spike_detector/Manual validation.xlsx'];

% Pt struct
data_folder = [locations.main_folder,'data/'];
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Get the indices of the patients with good spikes
T = readtable(validation_file);
good_pts = T.Var13;
good_pts = good_pts(~isnan(good_pts));
npts = length(good_pts);

count = 0;
% Loop over patients
for l = 1:npts
    j = good_pts(l);
    name = pt(j).name;
    
    %% Load the spike file
    fname = [spikes_folder,name,'_pc.mat'];
    if ~exist(fname,'file')
        fprintf('\nCannot find spike file for %s, skipping...\n',name);
        continue
    end
    
     %% Get basic info from the patient
    % load the spike file
    pc = load(fname);
    pc = pc.pc;
    
    % Skip the patient if it's incomplete
    if length(pc.file) < length(pt(j).ieeg.file) || ...
            length(pc.file(end).run) < size(pt(j).ieeg.file(end).run_times,1)
        fprintf('\n%s incomplete, skipping\n',name);
        continue
    end
        
    % add count
    count = count + 1;
    
    % reconcile files (deal with changes in electrode names)
    out = net_over_time(pc);
    out = reconcile_files(out);
    
    % Get the spikes and the labels
    times = out.times;
    spikes = out.montage(m).spikes;
    labels = out.montage(m).labels;
    
    % Clean the labels
    clean_labels = decompose_labels(labels,name);   
    
    % Get number of electrode localizations
    ne = length(pt(j).elecs);
    
    %% Find the anatomy corresponding to the spike labelsomy
    % Initialize a new spike_anatomy cell array corresponding to spike
    % labels
    spike_anatomy = cell(length(labels),1);
    spike_locs = nan(length(labels),3);
    already_filled = zeros(length(labels),1);
    
    for e = 1:ne
        % Get loc/anatomy names and labels
        ana_name = decompose_labels(pt(j).elecs(e).elec_names,name);
        locs = pt(j).elecs(e).locs;
        anatomy = pt(j).elecs(e).anatomy;

        % Indices of the loc/anatomy names that match the spike labels
        [lia,locb] = ismember(clean_labels,ana_name);
        % sanity check
        if ~isequal(clean_labels(lia~=0 & already_filled == 0),ana_name(locb(lia~=0 & already_filled == 0))), error('oh no'); end

        
        % Fill up spike anatomy and locs with the anatomy
        if ~strcmp(class(anatomy),'double')
            %spike_anatomy(lia~=0 & already_filled == 0) = anatomy(locb(lia~=0 & already_filled == 0));
            spike_anatomy(lia~=0) = anatomy(locb(lia~=0));
        end
        %spike_locs(lia~=0 & already_filled == 0,:) = locs(locb(lia~=0 & already_filled == 0),:);
        spike_locs(lia~=0,:) = locs(locb(lia~=0),:);
        
        % set already filled
        %already_filled(lia~=0 & already_filled == 0) = 1;
        already_filled(lia~=0) = 1;
    end
    
    % Get anatomical groupings
    [loc,lat] = cluster_anatomical_location(spike_anatomy);
    
    %% Put it all in the intermediate struct
    summ(count).name = name;
    summ(count).times = times;
    summ(count).spikes = spikes;
    summ(count).labels = clean_labels;
    summ(count).locs = spike_locs;
    summ(count).anatomy = spike_anatomy;
    summ(count).ana_loc = loc;
    summ(count).ana_lat = lat;


end

%% Save it all
save([out_folder,'summ.mat'],'summ');


