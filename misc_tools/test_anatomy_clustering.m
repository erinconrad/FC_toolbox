function test_anatomy_clustering(summ)

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Get the indices of the patients with good spikes
npts = length(summ);

% Loop over patients
for l = 1:npts

    % Get the spikes and the labels
    name = summ(l).name;
    labels = summ(l).labels;
    anatomy = summ(l).anatomy;
    loc = summ(l).ana_loc;
    lat = summ(l).ana_lat;
    
    % Show a table of the anatomy and assigned localization and
    % lateralization
    name
    table(labels,anatomy,lat,loc)
    pause
    
    
    
end


end