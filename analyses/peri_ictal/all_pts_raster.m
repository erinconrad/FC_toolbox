function all_pts_raster

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [results_folder,'peri_ictal_out/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Get patients
listing = dir([data_folder,'*.mat']);

for l = 1:length(listing)
    % load the file
    pc = load([data_folder,listing(l).name]);
    pc = pc.pc;
    
    % loop over files
    for f = 1:length(pc.file)
        % loop over seizures
        for s = 1:length(pc.file(f).sz)
            % do raster plot and save
            raster_peri_ictal(pc,f,s,1);
        end
        
    end
end

end