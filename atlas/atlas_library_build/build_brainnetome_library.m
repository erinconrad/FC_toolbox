function build_brainnetome_library

%% Parameters
atlas = 'brainnetome';

%% Get file locs
locations = fc_toolbox_locs;
atlas_folder = [locations.main_folder,'data/atlas/brainnetome/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


nums = [];
names = {};

%% Load the brainnetome list of regions
T =readtable([atlas_folder,'BNA_subregions.xlsx']);

lID = T.LabelID_L;
rID = T.LabelID_R;
reg_names = T.LeftAndRightHemisphere;

nrows = size(T,1);
for r = 1:nrows
    nums = [nums;lID(r);rID(r)]; % add the left and then the right
    names = [names;...
        [reg_names{r},'_L'];[reg_names{r},'_R']]; % add name for left and then right
    
end

%% Get unique ones
[nums,ia] = unique(nums);
names = names(ia);

%% Sort by number
[nums,I] = sort(nums);
names = names(I);

T = table(nums,names);
writetable(T,[atlas_folder,atlas,'_library.csv']);


end