function mt = find_mesial_temporal(names,which_atlas)

%% Get locs and lats for atlas names
names(cellfun(@isempty,names)) = {' '};
[~,lats] = lateralize_regions(names,which_atlas);
left = strcmp(lats,'L');
right = strcmp(lats,'R');

%% Define names corresponding to mesial temporal
switch which_atlas
    case 'aal'
        mt_names = {'Hippocampus','Amygdala'};
    case 'brainnetome'
        mt_names = {'Amyg','Hipp'};
end

mt = contains(names,mt_names);

%% Get mt L and R
mt_left = mt & left;
mt_right = mt & right;
mt_lr = [mt_left,mt_right];

end