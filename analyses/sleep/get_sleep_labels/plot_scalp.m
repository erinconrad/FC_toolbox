function plot_scalp(pt)

%% Parameters
% sample every how many blocks
block_stride = 24; % 4 hours = 24 ten-minute blocks

% Ignore this number of blocks at beginning and end of file
buffer = 36;% 6 hours = 36 ten-minute blocks

% scalp channels I want
scalp_labels = {'F3';'C3';'FZ';'CZ';'F4';'C4'};
montage = [1 2;3 4;5 6];

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Find pts with scalp
indices = find_pts_with_scalp(pt);

%% Loop over pts
for idx = 1:length(indices)
    p = indices(idx);
    name = pt(p).name;
    
    fprintf('\nDoing %s\n',name);
     
    % loop over files
    for f = 1:length(pt(p).ieeg.file)
        
        %% Check if I have the channels I want
        % DO THIS!!
        
        nblocks = size(pt(p).ieeg.file(f).block_times,1);
        
        %% loop over blocks (skipping buffer, striding)
        for b = buffer:block_stride:nblocks-buffer
            
            % Get run times
            run_times = run_times(b,:);
            
            
            
        end
    end
    
end


end