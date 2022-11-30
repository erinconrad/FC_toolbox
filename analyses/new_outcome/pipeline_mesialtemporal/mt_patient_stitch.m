function mt_patient_stitch(edf_path,name)


%% Load meta
meta = load([edf_path,'/meta.mat']);
meta = meta.meta;
times = meta.times;

%% Get allowable labels
allowable_labels = get_allowable_elecs;

%% Find labels that match allowable electrodes
allowed = ismember(meta.labels,allowable_labels);
allowed_labels = meta.labels(allowed);
nallowed = sum(allowed);

% Loop over files
nfiles = 72;
all_times = nan(nfiles,2);
all_bp = nan(nfiles,nallowed,5);
all_spikes = nan(nfiles,nallowed);
all_pc = nan(nfiles,nallowed,nallowed);
all_coh = nan(nfiles,nallowed,nallowed,6);
all_ad = nan(nfiles,nallowed);
all_rl = nan(nfiles,nallowed);
all_spike_times = cell(nfiles,1);

for f = 1:nfiles
    file_path = [edf_path,sprintf('/file%d.edf',f)];
    tic
    fprintf('\nDoing %s file %d of %d...',name,f,nfiles);

    %% Do the individual run
    out = individual_run_mt(file_path);
    fprintf('took %1.1f s\n',toc);

    %% Figure out times
    file_times = out.times;
    abs_times = times(f,1) + file_times;
    all_times(f,:) = abs_times;

    %% Figure out labels
    assert(equal(allowed_labels,out.clean_labels))
    
    %% Stitch together
    all_bp(f,:) = out.bp;
    spike_counts = nan(nallowed,1);
    %spike_counts = sum(out.spikes(:,1))

end

end