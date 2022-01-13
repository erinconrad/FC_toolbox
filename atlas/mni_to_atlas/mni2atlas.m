function mni_labels = mni2atlas(mni,atlas_path)

%% Read in atlas to get transformation matrix
V=niftiinfo(atlas_path); % get header
atlas = niftiread(V); % get 3D matrix
T=V.Transform.T; % get transformation matrix

%% Get coordinates of mni_locs in atlas image
mni_coords = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]*(inv(T));
mni_coords(:,4) = [];
mni_coords = round(mni_coords);

NN_flag = zeros(size(mni_coords,1), 2); % variable indicating no label initial found in atlas

%% Get atlas label based on image coordinates
mni_labels = nan(size(mni_coords,1),1);
for i = 1:size(mni_coords,1)

    if isnan(mni_coords(i,1)) || sum(mni_coords(i,:)>0) < 3 ||  mni_coords(i,1) > size(atlas,1) || ...
            mni_coords(i,2) > size(atlas,2) || mni_coords(i,3) > size(atlas,3)
        continue
    end
    
    mni_labels(i) = atlas(mni_coords(i,1), mni_coords(i,2), mni_coords(i,3));
    
    radius = 0; % initial radius size
    while mni_labels(i) == 0 % get mode of cubes until answer is achieved
        NN_flag(i,1) = 1; % store value indicating label not localizaed
        radius = radius + 1; % increase radius
        [my_sphere,coords] = gen_sphere(radius); % get coordinates to sample
        x = coords.x + mni_coords(i,1); 
        y = coords.y + mni_coords(i,2);
        z = coords.z + mni_coords(i,3);
        size(atlas);
        try
            ind = sub2ind(size(atlas),x,y,z); % convert to indices
            sphere_vals = atlas(ind); % sample sphere
            [~,~,vals] = find(sphere_vals); % find nonzero and non-NaN values
            if vals % if there are nonzero values, update label
                mni_labels(i) = mode(vals(:)); % get mode
                NN_flag(i,2) = radius; % store distance to NN
            end
        catch
            mni_labels(i) = 0; % get mode
            NN_flag(i,2) = radius; % store distance to NN
        end
        
        if radius>=6
            break
        end
        
    end
end

end

function [my_sphere,coords] = gen_sphere(radius)
    dim = 2*radius + 1; % size of box containing sphere
    mid = radius + 1; % midpoint of the box ([mid mid mid] is centroid)
    my_sphere = zeros(dim,dim,dim); % build box
    for i = 1:numel(my_sphere) % loop through each location in box
        [X,Y,Z] = ind2sub([dim,dim,dim],i); % get coordinates
            
        %% Erin switched this because I don't have "dist"
        %D = dist([mid mid mid; X Y Z]'); % distance from centroid
        D = vecnorm([mid mid mid] - [X Y Z]);
        
        if D <= radius % if less than radius, put in box
            my_sphere(i) = 1;
        end
    end
    
    % get coordinates and offset by center voxel
    [coords.x,coords.y,coords.z] = ind2sub(size(my_sphere), find(my_sphere));
    coords.x  = coords.x - mid;
    coords.y  = coords.y - mid;
    coords.z  = coords.z - mid;
end


