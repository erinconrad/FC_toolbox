function out = return_elec_locs(name,elec_loc_folder)

%% Parameters
%elec_file = 'electrodenames_coordinates_mni.csv';
elec_file = 'electrodenames_coordinates_native.csv';

listing = dir(elec_loc_folder);
match_idx = [];

for i = 1:length(listing)
    if contains(listing(i).name,name)
        match_idx = [match_idx;i];
    end
end

if length(match_idx) == 0
    error('no file match');
elseif length(match_idx) > 1
    error('multiple file matches');
end

%% Load the desired file
T = readtable([elec_loc_folder,listing(match_idx).name,'/',elec_file]);

%% Get electrode names
elec_names = T.Var1;
anatomy = T.Var2;
locs = [T.Var3 T.Var4 T.Var5];

%% Do some sanity checks
% Are the elec names strings
if ~iscell(elec_names)
    error('check elec names');
end

if ~strcmp(class(elec_names{1}),'char')
    error('check elec names');
end

% Are the anatomical names strings
if ~iscell(anatomy)
    error('check anatomy');
end

% Is the median distance between electrodes and the one listed under them
% close to 5 mm?
diff_locs = diff(locs,[],1);
dist_locs = vecnorm(diff_locs,2,2);
if abs(median(dist_locs)-5) > 0.1
    error('check distances');
end

out.elec_names = elec_names;
out.locs = locs;
out.anatomy = anatomy;

end