function out = reconcile_files(out)

%% Get all possible label names
nfiles = length(out.file);
all_labels = cell(nfiles,1);
n_runs_total = 0;
nmontages = length(out.file(1).montage);


for f = 1:nfiles
    labels = out.file(f).clean_labels;  
    all_labels{f} = labels;
    n_runs_total = n_runs_total + length(out.file(f).run_center);
end

%% Put all labels into one big list
all_labels_together = vertcat(all_labels{:});
unique_labels = unique(all_labels_together,'stable');
nchs = length(unique_labels);
out.all_labels = unique_labels;
out.n_runs_total = n_runs_total;
out.run_center = nan(n_runs_total,1);
out.file_index = nan(n_runs_total,1);
out.times = nan(n_runs_total,1);
out.mod_midnight = nan(n_runs_total,1);
nfreqs = size(out.file(1).montage(2).coh,2);

for m  = 1:nmontages
    %out.montage(m).net = nan(nchs,nchs,n_runs_total);
    out.montage(m).net = nan(nchs*(nchs-1)/2,n_runs_total);
    out.montage(m).spikes = nan(nchs,n_runs_total);
    out.montage(m).ad = nan(nchs,n_runs_total);
    %out.montage(m).coi_ch = nan(nchs,n_runs_total);
    out.montage(m).coi_global = nan(n_runs_total,1);
    out.montage(m).coa = nan(nchs*(nchs-1)/2,n_runs_total);
    out.montage(m).rl = nan(nchs,n_runs_total);
    out.montage(m).labels = cell(nchs,1);
    out.montage(m).n_rm_ictal = 0;
    out.montage(m).seq_info = nan(2,n_runs_total);
    out.montage(m).leader_montage = nan(nchs,n_runs_total);
    out.montage(m).coh = nan(nchs*(nchs-1)/2,nfreqs,n_runs_total);
    out.montage(m).bp = nan(nchs,nfreqs,n_runs_total);
end

%% Find how to match labels in each file with these labels
for f = 1:nfiles
    labels = out.file(f).clean_labels;  
    
    [lia,locb] = ismember(labels,unique_labels);
    
    % confirm that every member of labels is in unique labels
    if sum(lia) ~= length(labels), error('oh no'); end
    
    % locb contains index of unique_labels corresponding to each value of
    % labels - confirm!
    if ~isequal(unique_labels(locb),labels), error('oh no'); end
    
    out.file(f).full_idx = locb;
    
end

%% Go through and fill up data according to these indices
curr_run_idx = 1;
last_run_center = 0;
out.sz_times = [];
for f = 1:nfiles
    
    run_center = out.file(f).run_center;
    nruns = length(run_center);
    sz_times = out.file(f).sz_times;
    
    % file start time
    file_start_time = out.file(f).file_start_time;
    file_start_time = convert_prop_day_to_seconds_from_midnight(file_start_time);
    
    % get times of runs relative to seconds from midnight
    run_center_rel_midnight = run_center + file_start_time;
    
    % mod it
    run_center_rel_midnight = mod_secs_rel_24hrs(run_center_rel_midnight);
    out.mod_midnight(curr_run_idx:curr_run_idx+nruns-1,1) = run_center_rel_midnight;
    
    % run_center
    out.run_center(curr_run_idx:curr_run_idx+nruns-1,1) = run_center;
    out.file_index(curr_run_idx:curr_run_idx+nruns-1,1) = f;
    out.times(curr_run_idx:curr_run_idx+nruns-1) = run_center + last_run_center;
    out.sz_times = [out.sz_times;...
        sz_times + repmat(last_run_center,size(sz_times,1),1)];
    last_run_center = run_center(end)+last_run_center;
    
   
    % get indices in full set
    locb = out.file(f).full_idx;
    
    for m = 1:nmontages
        
        % get current data
        net = out.file(f).montage(m).net;
        spikes = out.file(f).montage(m).spikes;
        ad = out.file(f).montage(m).ad;
        coa = out.file(f).montage(m).coa;
        rl = out.file(f).montage(m).rl;
        %coi_ch = out.file(f).montage(m).coi_ch;
        coi_global = out.file(f).montage(m).coi_global;
        labels = out.file(f).montage(m).labels;
        n_rm_ictal = out.file(f).montage(m).n_rm_ictal;
        seq_info = out.file(f).montage(m).seq_info;
        leader_montage = out.file(f).montage(m).leader_montage;
        coh = out.file(f).montage(m).coh;
        bp = out.file(f).montage(m).bp;
        
        % prep net_uw and coa_uw
        net_uw = nan(length(labels),length(labels),nruns);
        coa_uw = nan(length(labels),length(labels),nruns);
        coh_uw = nan(length(labels),length(labels),nfreqs,nruns);
        
        % unwrap net
        for r = 1:nruns
            net_uw_temp = wrap_or_unwrap_adjacency_fc_toolbox(net(:,r));
            net_uw(:,:,r) = net_uw_temp;
            
            coa_uw_temp = wrap_or_unwrap_adjacency_fc_toolbox(coa(:,r));
            coa_uw(:,:,r) = coa_uw_temp;
            
            
            for i_f = 1:nfreqs
                coh_uw(:,:,i_f,r) = wrap_or_unwrap_adjacency_fc_toolbox(coh(:,i_f,r));
            end
            
        end
        
        % prep new ones (size of full set of labels)
        new_net_uw = nan(nchs,nchs,nruns);
        new_coa_uw = nan(nchs,nchs,nruns);
        new_coh_uw = nan(nchs,nchs,nfreqs,nruns);
        new_spikes = nan(nchs,nruns);
        new_ad = nan(nchs,nruns);
        new_bp = nan(nchs,nfreqs,nruns);
        %new_coi_ch = nan(nchs,nruns);
        new_rl = nan(nchs,nruns);
        new_labels = cell(nchs,1);
        new_leader_montage = nan(nchs,nruns);
        
        % fill the new ones based on the indices
        new_net_uw(locb,locb,:) = net_uw;
        new_coh_uw(locb,locb,:,:) = coh_uw;
        new_spikes(locb,:) = spikes;
        new_ad(locb,:) = ad;
        new_bp(locb,:,:) = bp;
        new_coa_uw(locb,locb,:) = coa_uw;
        new_rl(locb,:) = rl;
        new_leader_montage(locb,:) = leader_montage;
        %new_coi_ch(locb,:) = coi_ch;
        %new_ns(locb,:) = nansum(net_uw,2);
        
        % rewrap the net and coa
        new_net = nan(nchs*(nchs-1)/2,nruns);
        new_coa = nan(nchs*(nchs-1)/2,nruns);
        new_coh = nan(nchs*(nchs-1)/2,nfreqs,nruns);
        for r = 1:nruns
            new_net_temp = wrap_or_unwrap_adjacency_fc_toolbox(new_net_uw(:,:,r));
            new_net(:,r) = new_net_temp;
            
            
            new_coa_temp = wrap_or_unwrap_adjacency_fc_toolbox(new_coa_uw(:,:,r));
            new_coa(:,r) = new_coa_temp;
            
            for i_f = 1:nfreqs
                new_coh(:,i_f,r) = wrap_or_unwrap_adjacency_fc_toolbox(new_coh_uw(:,:,i_f,r));
            end
        end
        
        % fill the struct
        out.montage(m).net(:,curr_run_idx:curr_run_idx+nruns-1) = new_net;
        out.montage(m).coh(:,:,curr_run_idx:curr_run_idx+nruns-1) = new_coh;
        out.montage(m).spikes(:,curr_run_idx:curr_run_idx+nruns-1) = new_spikes;
        out.montage(m).ad(:,curr_run_idx:curr_run_idx+nruns-1) = new_ad;
        out.montage(m).bp(:,:,curr_run_idx:curr_run_idx+nruns-1) = new_bp;
        out.montage(m).coa(:,curr_run_idx:curr_run_idx+nruns-1) = new_coa;
        out.montage(m).rl(:,curr_run_idx:curr_run_idx+nruns-1) = new_rl;
        %out.montage(m).coi_ch(:,curr_run_idx:curr_run_idx+nruns-1) = new_coi_ch;
        out.montage(m).coi_global(curr_run_idx:curr_run_idx+nruns-1) = coi_global;
        out.montage(m).labels(locb) = labels;
        out.montage(m).n_rm_ictal = out.montage(m).n_rm_ictal + n_rm_ictal;
        out.montage(m).seq_info(:,curr_run_idx:curr_run_idx+nruns-1) = seq_info;
        out.montage(m).leader_montage(:,curr_run_idx:curr_run_idx+nruns-1) = new_leader_montage;
    end
    
    
    % go up on curr_run_idx
    curr_run_idx = curr_run_idx+nruns;
    
    
    
end

if 0
    im = 1;
    labels = out.montage(im).labels;
    all_nan = sum(~isnan(ns),2) == 0;
    turn_nans_gray(ns(~all_nan,:))
    yticks(1:sum(~all_nan));
    yticklabels(labels(~all_nan))
end


end