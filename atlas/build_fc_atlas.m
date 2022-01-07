function build_fc_atlas

%% To do
%{
- Add something to remove SOZ?
%}

%% Parameters
atlas = 'aal';
net_m = 2;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/atlas/'];
int_folder = [results_folder,'analysis/intermediate/'];
data_folder = [locations.main_folder,'data/'];
atlas_folder = [data_folder,'atlas/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load atlas library
library = readtable([atlas_folder,atlas,'_library.csv']);
atlas_nums = library.nums;
atlas_names = library.names;
n_parcels = length(atlas_nums);

%% Load pt folder
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

names = cell(npts,1);
sozs = cell(npts,1);

%% Initialize atlas
atlas_mat = nan(n_parcels,n_parcels,npts);

%% Loop over patients
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    name = summ.name;
    names{p} = name;
    
    %% Find corresponding pt index
    found_it = 0;
    for ip = 1:length(pt)
        if strcmp(pt(ip).name,name)
            found_it = 1;
            break
        end
    end
    if ~found_it, error('what'); end
    
    %% get rid and labels
    rid = pt(ip).rid;
    elabels = summ.labels;
    
    %% Get atlas
    out = get_atlas_parcellations(rid,elabels,name);
    enum = out.enum;
    
    %% Skip if empty
    if isempty(out.enum)
        continue
    end
    
    %% Get avg fc
    fc = summ.avg_fc;
    
    %% Soz
    soz = decompose_labels(summ.soz.labels,name);
    soz_out = get_atlas_parcellations(rid,soz,name);
    soz_num = soz_out.enum;
    sozs{p} = soz_num;
    
    %% Put into atlas space
    fc_atlas_space = nan(n_parcels,n_parcels);
    
    % Loop over n_parcels
    for inp = 1:n_parcels
        
        % get the atlas identifier
        which_enum_i = atlas_nums(inp);
        
        % which elecs for this patient belong to that
        which_elecs_i = enum == which_enum_i;
        
        % Loop over again
        for jnp = 1:inp-1
            
            % get the atlas identifier
            which_enum_j = atlas_nums(jnp);

            % which elecs for this patient belong to that
            which_elecs_j = enum == which_enum_j;
            
            % Assign the average of all the functional connectivities
            % matching these to be the fc edge
            fc_atlas_space(inp,jnp) = nanmean(fc(which_elecs_i,which_elecs_j),'all');
            fc_atlas_space(jnp,inp) = nanmean(fc(which_elecs_i,which_elecs_j),'all');
            
        end
        
    end
    
    atlas_mat(:,:,p) = fc_atlas_space;
    
end

%% Save the atlas
out.atlas = atlas_mat;
out.atlas_nums = atlas_nums;
out.atlas_names = atlas_names;
out.pt_names = names;
out.sozs = sozs;
save([out_folder,atlas,'.mat'],'out');

end