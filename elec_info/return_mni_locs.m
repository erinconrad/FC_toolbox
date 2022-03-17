function out = return_mni_locs(name,elec_loc_folder)

%% Parameters
elec_file = 'electrodenames_coordinates_mni.csv';
%elec_file = 'electrodenames_coordinates_native.csv';

listing = dir(elec_loc_folder);
match_idx = [];

for i = 1:length(listing)
    if contains(listing(i).name,name)
        match_idx = [match_idx;i];
    end
end

if isempty(match_idx)
    fprintf('\nWarning, no folder match for %s\n',name);
    out = [];
    return
else
    fprintf('\nGetting elec locs for %s\n',name);
end

for i = 1:length(match_idx)
    
    which_index = match_idx(i);

    %% Load the desired file
    try
        T = readtable([elec_loc_folder,listing(which_index).name,'/',elec_file]);
    catch
        try
            
            folder_name = listing(which_index).name;
            C = strsplit(folder_name,'_');
            rid = C{1};
            T = readtable([elec_loc_folder,listing(which_index).name,'/',rid,'/',elec_file]);
        catch
            fprintf('\nWarning, no file match for %s\n',name);

            out(i).folder_name = listing(which_index).name;
            out(i).elec_names = [];
            out(i).locs = [];
            out(i).anatomy = [];

            continue
        end
    end
    
    if ~ismember(T.Properties.VariableNames,'Var1')
        out(i).folder_name = listing(which_index).name;
        out(i).elec_names = [];
        out(i).locs = [];
        out(i).anatomy = [];
        continue
    end

    %% Get electrode names
    elec_names = T.Var1;
    %anatomy = T.Var2;
    locs = [T.Var2 T.Var3 T.Var4];

    %% Do some sanity checks
    % Are the elec names strings
    if ~iscell(elec_names)
        error('check elec names');
    end

    if ~strcmp(class(elec_names{1}),'char')
        error('check elec names');
    end
    
    % Is the median distance between electrodes and the one listed under them
    % close to 5 mm?
    %{
    diff_locs = diff(locs,[],1);
    dist_locs = vecnorm(diff_locs,2,2);
    if abs(median(dist_locs)-5) > 0.1
        error('check distances');
    end
    %}
    out(i).folder_name = listing(which_index).name;
    out(i).elec_names = elec_names;
    out(i).locs = locs;
    %out(i).anatomy = anatomy;
end

end