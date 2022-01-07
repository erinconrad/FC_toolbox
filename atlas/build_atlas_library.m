function build_atlas_library

%% Parameters
atlas = 'mni';

%% Get file locs
locations = fc_toolbox_locs;
atlas_folder = [locations.main_folder,'data/atlas/'];
parcellation_folder = [atlas_folder,'atlas_localizations/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

listing = dir([parcellation_folder,'*.csv']);

nums = [];
names = {};

%% Loop over all atlas files
for l = 1:length(listing)
    T = readtable([parcellation_folder,listing(l).name]);
    
    %% Choose the atlas I want
    switch atlas
        case 'aal'
            anum = T.AAL_region_number;
            aname = T.AAL_label;
        case 'mni'
            anum = T.mni_icbm152_CerebrA_tal_nlin_sym_09c_region_number;
            aname = T.mni_icbm152_CerebrA_tal_nlin_sym_09c_label;
    end
    
    %% Find those that aren't in the library already
    lia = ~ismember(anum,nums);
    
    %% Add those that aren't in the library
    nums = [nums;anum(lia)];
    names = [names;aname(lia)];
    
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