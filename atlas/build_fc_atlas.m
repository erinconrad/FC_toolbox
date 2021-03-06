function build_fc_atlas

%% Parameters
atlas = 'aal_bernabei';%'brainnetome'; %'aal_bernabei'; %
too_many_spikes = 1; % 1 spikes/elecs/min

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
n_elecs_all = nan(n_parcels,npts);
elecs_atlas = cell(npts,1);
elecs_labels = cell(npts,1);
spikes_atlas = nan(n_parcels,npts);
spike_leader_atlas = nan(n_parcels,npts);
normal_atlas_mat = nan(n_parcels,n_parcels,npts);
normal_include = nan(n_parcels,npts);
all_soz_locs = cell(npts,1);
all_soz_lats = cell(npts,1);

%% Loop over patients
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    name = summ.name;
    names{p} = name;
    good_spikes = summ.good_spikes;
    
    %% Find corresponding pt index
    found_it = 0;
    for ip = 1:length(pt)
        if strcmp(pt(ip).name,name)
            found_it = 1;
            break
        end
    end
    if ~found_it, error('what'); end
    
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
        fc = summ.avg_fc;
    end
    
    elecs_labels{p} = elabels;
    
    %% Get spikes
    if good_spikes
        spikes = summ.spikes;
        leader = summ.leader;
    else
        spikes = nan(size(summ.spikes));
        leader = nan(size(summ.leader));
    end
    
    %% Find and remove non intracranial
    ekg = find_non_intracranial(elabels);
    fc(ekg,:) = nan;
    fc(:,ekg) = nan;
    spikes(ekg,:) = nan;
    leader(ekg,:) = nan;
    
    %% Get average spikes over time
    spikes = nanmean(spikes,2);
    leader = nanmean(leader,2);
    

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
    soz_loc = summ.soz.loc;
    soz_lat = summ.soz.lat;
    all_soz_locs{p} = soz_loc;
    all_soz_lats{p} = soz_lat;
    
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
    
    %% Find channels that I want to exclude from the normal atlas
    % Specifically, SOZ channels and channels with too many spikes
    nchs = length(spikes);
    is_soz = zeros(nchs,1);
    is_soz(soz) = 1;
    exclude = is_soz | spikes > too_many_spikes;
    
    soz_num = out.enum(soz);
    sozs{p} = soz_num;
    
    %% Put into atlas space
    fc_atlas_space = nan(n_parcels,n_parcels);
    n_elecs = nan(n_parcels,1);
    electrode_atlas_assignment = nan(size(fc,1),1);
    norm_fc_atlas_space = nan(n_parcels,n_parcels);
    
    % Loop over n_parcels
    for inp = 1:n_parcels
        
        % get the atlas identifier
        which_enum_i = atlas_nums(inp);
        
        % which elecs for this patient belong to that
        which_elecs_i = out.enum == which_enum_i;
        
        % Get number of electrodes contributing to parcel
        n_elecs(inp) = sum(which_elecs_i);
        
        % Get average spikes and leader
        spikes_atlas(inp,p) = nanmean(spikes(which_elecs_i));
        spike_leader_atlas(inp,p) = nanmean(leader(which_elecs_i));
        
        normal_include(inp,p) = sum(which_elecs_i & ~exclude);
        
        % Loop over again
        for jnp = 1:inp-1
            
            % get the atlas identifier
            which_enum_j = atlas_nums(jnp);

            % which elecs for this patient belong to that
            which_elecs_j = out.enum == which_enum_j;
            electrode_atlas_assignment(which_elecs_j) = which_enum_j;
            
            % Assign the average of all the functional connectivities
            % matching these to be the fc edge
            fc_atlas_space(inp,jnp) = nanmean(fc(which_elecs_i,which_elecs_j),'all');
            fc_atlas_space(jnp,inp) = nanmean(fc(which_elecs_i,which_elecs_j),'all');
            
            % Now, make a normal atlas where I exclude bad channels
            norm_fc_atlas_space(inp,jnp) = nanmean(fc(which_elecs_i & ~exclude,which_elecs_j & ~ exclude),'all');
            norm_fc_atlas_space(jnp,inp) = nanmean(fc(which_elecs_i & ~exclude,which_elecs_j & ~ exclude),'all');
            
            
            
        end
        
    end
    
    atlas_mat(:,:,p) = fc_atlas_space;
    n_elecs_all(:,p) = n_elecs;
    elecs_atlas{p} = electrode_atlas_assignment;
    normal_atlas_mat(:,:,p) = norm_fc_atlas_space;
    
end

%% Save the atlas
out.atlas = atlas_mat;
out.normal_atlas = normal_atlas_mat;
out.atlas_nums = atlas_nums;
out.atlas_names = atlas_names;
out.n_elecs_all = n_elecs_all;
out.normal_include = normal_include;
out.pt_names = names;
out.sozs = sozs;
out.missing_names = missing_names;
out.elecs_atlas = elecs_atlas;
out.elecs_labels = elecs_labels;
out.spikes_atlas = spikes_atlas;
out.spike_leader_atlas = spike_leader_atlas;
out.all_soz_locs = all_soz_locs;
out.all_soz_lats = all_soz_lats;
save([out_folder,atlas,'.mat'],'out');


end