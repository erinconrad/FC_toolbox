function plot_all_sleep_rasters

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
raster_out_folder = [out_folder,'rasters/'];
if ~exist(raster_out_folder,'dir')
    mkdir(raster_out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the out file
load([out_folder,'out.mat']);

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

for p = 2%1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    name = summ.name;
    sleep_raster(summ,out,p);
    print(gcf,[raster_out_folder,name],'-dpng')
    close all
    
end



end