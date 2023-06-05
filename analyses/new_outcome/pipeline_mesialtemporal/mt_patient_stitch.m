function mt_patient_stitch(pt,which_pt,edf_path,edf_summ_path,name,overwrite,overlap_log_file,szT,mT,attempt_remove_oob)

%% Seed random number generator so that I get the same result each time I do this.
rng(0)


%% See if I've done it
if exist([edf_summ_path,name,'/summ.mat'],'file')~=0
    if overwrite == 0
        fprintf('\nAlready did %s, skipping.\n',name);
        return
    else
        fprintf('\nOverwriting %s.\n',name);
    end
end

%% Load meta
if ~exist([edf_path,name,'/meta.mat'],'file')
    return
end
meta = load([edf_path,name,'/meta.mat']);
meta = meta.meta;
times = meta.times;

%% Find labels that match allowable electrodes
% Load first file to get labels
info = edfinfo([edf_path,name,sprintf('/file%d.edf',1)]);
labels = cellstr(info.SignalLabels);
%allowed = ismember(labels,allowable_labels);
%allowed_labels = labels(allowed);
%nallowed = sum(allowed);

% Get allowable labels
name = pt(which_pt).name;
potentially_allowable_labels = get_allowable_elecs(name);

%% Find which of these electrodes are outside the brain
if isempty(pt(which_pt).dkt)

else
    % Get atlas labels
    dkt = pt(which_pt).dkt.label;
    atropos = pt(which_pt).atropos.label;
    atlas_elec_names = pt(which_pt).atropos.names;
    
    % convert atlas names to the A, B, C convention
    atlas_elec_names = mt_name_conversion(atlas_elec_names,name);
    
    if attempt_remove_oob
        outside_brain = zeros(length(potentially_allowable_labels),1);
        for i = 1:length(potentially_allowable_labels)
            % find the matching atlas elec name
            match = strcmp(potentially_allowable_labels{i},atlas_elec_names);
        
            if match == 0, continue; end
        
            if strcmp(dkt{match},'EmptyLabel') && (strcmp(atropos{match},'CSF') ...
                    || strcmp(atropos{match},'EmptyLabel'))
                outside_brain(i) = 1;
            end
        end
        outside_brain = logical(outside_brain);
        
        if 0
            table(pt(which_pt).atropos.names,atlas_elec_names,atropos,dkt)
        end
        
        % remove those outside brain
        potentially_allowable_labels(outside_brain) = [];
    end
end

% Remove some potentially allowable labels if they aren't really targeting mesial temporal region for that patient
exr = strcmp(mT.name,name); assert(sum(exr==1));
exc = mT.exclude{exr};
if ~isempty(exc)
    C = strsplit(exc,', ');
    rm_allow = zeros(length(potentially_allowable_labels),1);
    for i = 1:length(C) % loop over exclusion labels
        rm_allow(contains(potentially_allowable_labels,C{i},'ignorecase',true)) = 1;
    end
    potentially_allowable_labels(rm_allow==1) = [];
end

% Turn non-A, B, C MT electrodes into A, B, C
old_labels = labels;
labels = mt_name_conversion(labels,name);

% Find labels that match allowable electrodes and have symmetric coverage
allowed_labels = find_mt_symmetric_coverage(labels,potentially_allowable_labels);
if isempty(allowed_labels)
    fprintf('\n No allowed electrodes for %s, skipping.\n',name);
    return
end
nallowed = length(allowed_labels);

% Get demographic data
clinical = pt(which_pt).clinical;

if isempty(clinical)
    surgery = [];
    resec_lat = [];
    resec_loc = [];
    ablate_lat = [];
    ablate_loc = [];
    engel = [];
    ilae = [];

else
    surgery = clinical.surgery;
    resec_lat = clinical.resection_lat;
    resec_loc = clinical.resection_loc;
    ablate_lat = clinical.ablation_lat;
    ablate_loc = clinical.ablation_loc;
    engel = {clinical.engel{1},clinical.engel{2}};
    ilae = {clinical.ilae{1},clinical.ilae{2}};
end

% Get the correct row of the SOZ table
szr = nan;
for k = 1:size(szT,1)
    if strcmp(szT.name{k},name)
        szr = k;
        break
    end
end
soz_loc = szT.region{k};
soz_lat = szT.lateralization{k};


% Loop over files
nfiles = 72;
nmontages = 3;
all_times = nan(nfiles,2);
all_rel_bp = nan(nfiles,nmontages,nallowed,6);
all_bp = nan(nfiles,nmontages,nallowed,6);
all_spike_counts = nan(nfiles,nmontages,nallowed);
all_pc = nan(nfiles,nmontages,nallowed,nallowed);
all_pc_squared = nan(nfiles,nmontages,nallowed,nallowed);
all_coh = nan(nfiles,nmontages,nallowed,nallowed,6);
all_plv = nan(nfiles,nmontages,nallowed,nallowed,6);
all_re = nan(nfiles,nmontages,nallowed,nallowed,6);
all_ad = nan(nfiles,nmontages,nallowed);
all_rl = nan(nfiles,nmontages,nallowed);
all_se = nan(nfiles,nmontages,nallowed);
all_xcor = nan(nfiles,nmontages,nallowed,nallowed);
all_lags = nan(nfiles,nmontages,nallowed,nallowed);
all_spike_times = cell(nmontages,1);
all_ll = nan(nfiles,nmontages,nallowed);
montages = cell(nmontages,3);
montage_labels = cell(nmontages,3);
skipped_file = zeros(nfiles,1);

all_is_run = nan(nfiles,nmontages,nallowed);

skip_pt = 0;
for f = 1:nfiles
    file_path = [edf_path,name,sprintf('/file%d.edf',f)];
    tic
    fprintf('\nDoing %s file %d of %d...',name,f,nfiles);

    %% Do the individual run
    out = individual_run_mt(file_path,pt,which_pt,meta,f,overlap_log_file,mT,attempt_remove_oob);
    fprintf('took %1.1f s\n',toc);
    if isempty(out) 
        skipped_file(f) = 1;
        continue; 
    end


    %% Figure out times
    file_times = out.times;
    abs_times = times(f,1) + file_times;
    all_times(f,:) = abs_times;

    %% Figure out labels
    % Convert MUSC labels into HUP labels for consistency
    allowed_labels = convert_musc_labels_to_hup(allowed_labels);
    assert(isequal(allowed_labels,out.clean_labels))
    
    for im = 1:nmontages
        %% Stitch together
        % Get spikes
        gdf = out.montage(im).gdf;
        
        if ~isempty(gdf)
        % re-align index to file index
            gdf(:,2) = gdf(:,2) + out.idx(1) -1; % if gdf index is 1, that means it occurs at rand_start;
        end
    
        % save spike times (for future error checking and plotting)
        all_spike_times{im} = [all_spike_times{im};gdf repmat(f,size(gdf,1),1)];
    
        % get spike counts
        if ~isempty(gdf)
            X = gdf(:,1); % get channels
            spike_counts = accumarray(X, ones(size(X)), [nallowed 1], @sum); % this gets the count of spikes for each channe;
            
            % make skipped channels nans rather than zeros
            spike_counts(out.montage(im).skip) = nan;
            all_spike_counts(f,im,:) = spike_counts;
        else
            all_spike_counts(f,im,:) = 0;
        end
    
        % get RL
        if ~isempty(gdf)
            timing = gdf(:,3);
            % take the mean timing for each channel
            rl = accumarray(X, timing, [nallowed 1], @mean,nan); % if no spikes, make this nan
            all_rl(f,im,:) = rl;
        else
            all_rl(f,im,:) = nan;
        end
    
        % get the other stuff
        all_bp(f,im,:,:) = out.montage(im).bp;
        all_rel_bp(f,im,:,:) = out.montage(im).rel_bp;
        all_pc(f,im,:,:) = out.montage(im).pc;
        all_pc_squared(f,im,:,:) = out.montage(im).pc_squared;
        all_se(f,im,:) = out.montage(im).se;
        all_xcor(f,im,:,:) = out.montage(im).xcor;
        all_lags(f,im,:,:) = out.montage(im).lags;
        all_coh(f,im,:,:,:) = out.montage(im).coh;
        all_plv(f,im,:,:,:) = out.montage(im).plv;
        all_re(f,im,:,:,:) = out.montage(im).re;
        all_ad(f,im,:) = out.montage(im).ad;
        all_is_run(f,im,:) = out.montage(im).is_run;
        all_ll(f,im,:) = out.montage(im).ll;

        montages{im} = out.montage(im).name;
        montage_labels{im} = out.montage(im).labels;
    end


end


%% Plot random spike detections
if ~exist([edf_summ_path,name],'dir')
    mkdir([edf_summ_path,name])
end


%% Output the stuff
nout.all_times = all_times;
nout.all_bp = all_bp;
nout.all_rel_bp = all_rel_bp;
nout.all_spike_counts = all_spike_counts;
nout.all_pc = all_pc;
nout.all_pc_squared = all_pc_squared;
nout.all_coh = all_coh;
nout.all_re = all_re;
nout.all_ad = all_ad;
nout.all_plv = all_plv;
nout.all_rl = all_rl;
nout.all_ll = all_ll;
nout.all_se = all_se;
nout.all_xcor = all_xcor;
nout.all_lags = all_lags;
nout.all_spike_times = all_spike_times;
nout.fs = out.fs;
nout.edf_path = edf_path;
nout.name = name;
nout.labels = out.clean_labels;
nout.montage_labels = montage_labels;
nout.all_is_run = all_is_run;
nout.montages = montages;
nout.surgery = surgery;
nout.resec_lat = resec_lat;
nout.resec_loc = resec_loc;
nout.ablate_lat = ablate_lat;
nout.ablate_loc = ablate_loc;
nout.engel = engel;
nout.ilae = ilae;
nout.soz_loc = soz_loc;
nout.soz_lat = soz_lat;
nout.allowed_labels = allowed_labels;
nout.old_allowed_labels = out.old_allowed_labels;
nout.skipped_file = skipped_file;

out = nout;

%% Save the file
if ~exist([edf_summ_path,name],'dir')
    mkdir([edf_summ_path,name])
end
save([edf_summ_path,name,'/summ.mat'],'out');

for im = 1:nmontages
    plot_random_spikes(all_spike_times{im},name,nout.labels,montages{im},edf_path,edf_summ_path,mT,pt,which_pt,attempt_remove_oob)
end

end