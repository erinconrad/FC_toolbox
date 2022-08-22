function mt_lr = return_bilateral_mesial_temporal_spikes(spikes,names,which_atlas)


%% Get locs and lats for atlas names
names(cellfun(@isempty,names)) = {' '};
[locs,lats] = lateralize_regions(names,which_atlas);
left = strcmp(lats,'L');
right = strcmp(lats,'R');
neither_lat = ~left & ~right;


%% Define names corresponding to mesial temporal
switch which_atlas
    case 'aal'
        mt_names = {'Hippocampus','Amygdala'};
    case 'brainnetome'
        mt_names = {'Amyg','Hipp'};
end

mt = contains(names,mt_names);

if 0
names(mt & left)
names(mt & right)
end

%% Get mt L and R
mt_left = nanmean(spikes(mt & left));
mt_right = nanmean(spikes(mt & right));
mt_lr = [mt_left,mt_right];

end