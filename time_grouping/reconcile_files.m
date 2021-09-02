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

for m  = 1:nmontages
    %out.montage(m).net = nan(nchs,nchs,n_runs_total);
    out.montage(m).net = nan(nchs*(nchs-1)/2,n_runs_total);
    out.montage(m).spikes = nan(nchs,n_runs_total);
    out.montage(m).ad = nan(nchs,n_runs_total);
    out.montage(m).labels = cell(nchs,1);
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
for f = 1:nfiles
    
    run_center = out.file(f).run_center;
    nruns = length(run_center);
    
    % run_center
    out.run_center(curr_run_idx:curr_run_idx+nruns-1,1) = run_center;
    out.file_index(curr_run_idx:curr_run_idx+nruns-1,1) = f;
    out.times(curr_run_idx:curr_run_idx+nruns-1) = run_center + last_run_center;
    last_run_center = run_center(end)+last_run_center;
   
    % get indices in full set
    locb = out.file(f).full_idx;
    
    for m = 1:nmontages
        
        % get current data
        net = out.file(f).montage(m).net;
        spikes = out.file(f).montage(m).spikes;
        ad = out.file(f).montage(m).ad;
        labels = out.file(f).montage(m).labels;
        
        % prep net_uw
        net_uw = nan(length(labels),length(labels),nruns);
        
        % unwrap net
        for r = 1:nruns
            net_uw_temp = wrap_or_unwrap_adjacency_fc_toolbox(net(:,r));
            net_uw(:,:,r) = net_uw_temp;
        end
        
        % prep new ones (size of full set of labels)
        new_net_uw = nan(nchs,nchs,nruns);
        new_spikes = nan(nchs,nruns);
        new_ad = nan(nchs,nruns);
        new_labels = cell(nchs,1);
        
        % fill the new ones based on the indices
        new_net_uw(locb,locb,:) = net_uw;
        new_spikes(locb,:) = spikes;
        new_ad(locb,:) = ad;
        %new_ns(locb,:) = nansum(net_uw,2);
        
        % rewrap the net
        new_net = nan(nchs*(nchs-1)/2,nruns);
        for r = 1:nruns
            new_net_temp = wrap_or_unwrap_adjacency_fc_toolbox(new_net_uw(:,:,r));
            new_net(:,r) = new_net_temp;
        end
        
        % fill the struct
        out.montage(m).net(:,curr_run_idx:curr_run_idx+nruns-1) = new_net;
        out.montage(m).spikes(:,curr_run_idx:curr_run_idx+nruns-1) = new_spikes;
        out.montage(m).ad(:,curr_run_idx:curr_run_idx+nruns-1) = new_ad;
        out.montage(m).labels(locb) = labels;
        
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