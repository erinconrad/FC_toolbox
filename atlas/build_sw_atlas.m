function build_sw_atlas


%% Parameters
atlas = 'aal_bernabei';

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

sw_disc_folder = [scripts_folder,'analyses/sleep/data/'];

%% Get normalized ADR value that best discriminates sleep from wake
sw_out = load([sw_disc_folder,'out.mat']);
sw_out = sw_out.out;
sw_disc = sw_out.roc_out.disc;

%% Load atlas library
if contains(atlas,'aal_bernabei')
    [atlas_names,atlas_nums] = aal_region_to_name([atlas_folder,'AAL116_WM.txt'],[]);
elseif contains(atlas,'brainnetome')
    library = readtable([atlas_folder,atlas,'_library.xlsx']);
    atlas_nums = library.nums;
    atlas_names = library.names;
else
    library = readtable([atlas_folder,atlas,'_library.csv'],'ReadVariableNames',false);
    atlas_nums = library.nums;
    atlas_names = library.names;
    
end
n_parcels = length(atlas_nums);

%% Load pt folder
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

names = cell(npts,1);
sozs = cell(npts,1);

%% Missing pts
missing_names = {};

%% Initialize atlas
atlas_mat = nan(n_parcels,n_parcels,npts); 
atlas_mat_ws = nan(n_parcels,n_parcels,2,npts);  % wake and sleep
n_elecs_all = nan(n_parcels,npts);

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
    
    %% Get sleep and wake indices
    ad = summ.ad;
    labels = summ.labels;
    
    % remove intracranial
    ekg = find_non_intracranial(labels);
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    
    % Get sleep and wake times
    [sleep,wake] = find_sleep_wake(ad,[],sw_disc);
    
    %% get stuff
    rid = pt(ip).rid;
    
    if contains(atlas,'bipolar')
        elabels = summ.bipolar_labels;
        locs = summ.bipolar_locs;
        fc = summ.avg_fc_bi;
        bipolar_pair = summ.bipolar_pair;
    else
        elabels = summ.labels;
        locs = summ.locs;
        fc = summ.avg_fc_car;
        fc_all = summ.fc_car;
    end
    
        
    %% Get atlas
    if contains(atlas,'aal_bernabei')
        
        % John's method (convert loc in MNI space to AAL atlas) 
        
        [~, mni_roi, ~] = nifti_values_erin(locs,[atlas_folder,'AAL116_WM.nii']);
        %aal_names = aal_region_to_name([atlas_folder,'AAL116_WM.txt'],mni_roi);
        out.enum = mni_roi;
        out.ename = atlas_names;
    
    elseif contains(atlas,'brainnetome')
        mni_roi = mni2atlas(locs,[atlas_folder,'brainnetome/BN_Atlas_246_2mm.nii.gz']);
        out.enum = mni_roi;
        out.ename = atlas_names;
    else
        % Andy method (table with parcellations for each electrode)
        out = get_atlas_parcellations(rid,elabels,name);
      
    end
    
    if 0
        table(elabels,out.ename)
    end
        
    
    %% Skip if empty
    if isempty(out.enum)
        missing_names = [missing_names;name];
        continue
    end
    

    %% Soz
    %{
    soz = decompose_labels(summ.soz.labels,name);
    soz_out = get_atlas_parcellations(rid,soz,name);
    soz_num = soz_out.enum;
    sozs{p} = soz_num;
    %}
    soz = summ.soz.chs;
    soz(soz==0) = [];
    
    if contains(atlas,'bipolar')
        soz_bipolar = [];
        % find bipolar channels that contain these contacts
        for ib = 1:size(bipolar_pair,1)
            if any(ismember(soz,bipolar_pair(ib,:)))
                soz_bipolar = [soz_bipolar;ib];
            end
        end
        soz_bipolar = unique(soz_bipolar);
        soz = soz_bipolar;
    end
    
    
    soz_num = out.enum(soz);
    sozs{p} = soz_num;
    
    %% Put into atlas space
    fc_atlas_space = nan(n_parcels,n_parcels);
    fc_atlas_space_ws = nan(n_parcels,n_parcels,2); % wake and sleep
    n_elecs = nan(n_parcels,1);
    
    % Loop over n_parcels
    for inp = 1:n_parcels
        
        % get the atlas identifier
        which_enum_i = atlas_nums(inp);
        
        % which elecs for this patient belong to that
        which_elecs_i = out.enum == which_enum_i;
        
        % Get number of electrodes contributing to parcel
        n_elecs(inp) = sum(which_elecs_i);
        
        % Loop over again
        for jnp = 1:inp-1
            
            % get the atlas identifier
            which_enum_j = atlas_nums(jnp);

            % which elecs for this patient belong to that
            which_elecs_j = out.enum == which_enum_j;
            
            % Assign the average of all the functional connectivities
            % matching these to be the fc edge
            fc_atlas_space(inp,jnp) = nanmean(fc(which_elecs_i,which_elecs_j),'all');
            fc_atlas_space(jnp,inp) = nanmean(fc(which_elecs_j,which_elecs_i),'all');
            
            % Sleep wake
            fc_atlas_space_ws(inp,jnp,1) = nanmean(fc_all(which_elecs_i,which_elecs_j,wake),'all');
            fc_atlas_space_ws(jnp,inp,1) = nanmean(fc_all(which_elecs_j,which_elecs_i,wake),'all');
            fc_atlas_space_ws(inp,jnp,2) = nanmean(fc_all(which_elecs_i,which_elecs_j,sleep),'all');
            fc_atlas_space_ws(jnp,inp,2) = nanmean(fc_all(which_elecs_j,which_elecs_i,sleep),'all');
            
            
        end
        
    end
    
    atlas_mat(:,:,p) = fc_atlas_space;
    atlas_mat_ws(:,:,:,p) = fc_atlas_space_ws;
    n_elecs_all(:,p) = n_elecs;
    
end

%% Save the atlas
out.atlas = atlas_mat;
out.atlas_ws = atlas_mat_ws;
out.atlas_nums = atlas_nums;
out.atlas_names = atlas_names;
out.n_elecs_all = n_elecs_all;
out.pt_names = names;
out.sozs = sozs;
out.missing_names = missing_names;
save([out_folder,atlas,'_ws.mat'],'out');


end