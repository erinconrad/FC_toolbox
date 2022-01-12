function assign_aal_to_canon

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/atlas/'];
int_folder = [results_folder,'analysis/intermediate/'];
data_folder = [locations.main_folder,'data/'];
atlas_folder = [data_folder,'atlas/'];
yeo_canon_networks = [data_folder,'canonical_networks/Yeo_JNeurophysiol11_MNI152/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load aal image
V_aal=niftiinfo([atlas_folder,'AAL116_WM.nii']); % get header
atlas_aal = niftiread(V_aal); % get 3D matrix
T_aal=V_aal.Transform.T; % get transformation matrix
T_aal=T_aal'; % transpose transformation matrix

%% Load Yeo networks (liberal mask)
V_yeo=niftiinfo([yeo_canon_networks,'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz']); % get header
atlas_yeo = niftiread(V_yeo); % get 3D matrix
T_yeo=V_yeo.Transform.T; % get transformation matrix
T_yeo=T_yeo'; % transpose transformation matrix

%% Get the coordinates of each aal parcel
% convert size of aal atlas to an Nx3 matrix where N is the number of
% parcels
nparcels = size(atlas_aal,1)*size(atlas_aal,2)*size(atlas_aal,3);
parcel_coords = nan(nparcels,3);
aal_assignments = nan(nparcels,1);
count = 0;
for i = 1:size(atlas_aal,1)
    for j = 1:size(atlas_aal,2)
        for k = 1:size(atlas_aal,3)
            count = count + 1;
            parcel_coords(count,:) = [i j k];
            aal_assignments(count,:) = atlas_aal(i,j,k);
        end
    end
end
assert(count == nparcels)

% convert coords to mni
mni = cor2mni(parcel_coords, T_aal);

%% Get Yeo assignments of each mni
yeo_coordinates = mni2cor(mni, T_yeo);
yeo_assignments = nan(size(yeo_coordinates,1),1);
for i = 1:size(yeo_coordinates,1)
    yeo_assignments(i) = atlas_yeo(yeo_coordinates(i,1),yeo_coordinates(i,2),yeo_coordinates(i,3));
end
unique_yeo = unique(yeo_assignments);
n_yeo = length(unique_yeo);

%% Get all yeo assignments corresponding to each aal assignenment
%
unique_aal = unique(aal_assignments);
n_aal = length(unique_aal);
aal_to_yeo = nan(n_aal,n_yeo);
for i = 1:n_aal
    
    % find the indices belonging to that
    idx = aal_assignments == unique_aal(i);
    
    % What percent the coordinates belonging to this aal belong to each
    % yeo?
    for j = 1:n_yeo
        curr_yeo = unique_yeo(j);
        aal_to_yeo(i,j) = sum(yeo_assignments(idx) == curr_yeo)/sum(idx);
    end
    
end
%}


end