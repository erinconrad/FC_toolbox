function assignments = parcellate_atlas(elecs)

%% Parameters
dist_thresh = 1; % how close does the electrode need to be to be in the atlas

%% Load the atlas image
img = atlas.img;


%% Get electrode coordinates into whatever coordinate system we need for the atlas

%% Get coordinates of points in atlas
[xind,yind,zind] = ind2sub(size(img),find(ismember(img,1:max(max(max(img))))));
atlas_coordinates = [xind,yind,zind];

%% Assign electrodes to ROIs in array space
% Get closest atlas point for each electrode
[k,dist] = dsearchn(atlas_coordinates,elec_coordinates);

% which electrodes are close enough to an atlas point
close_enough = dist <= dist_thresh;

% Get the image 
assignments = nan(nelecs,1);
assignments(close_enough) = img(k(close_enough));

%% Show electrodes and assignments
if 1
    vertical_slice = 100;
    turn_nans_white(img(:,:,vertical_slice))
    hold on
    
    
end

end