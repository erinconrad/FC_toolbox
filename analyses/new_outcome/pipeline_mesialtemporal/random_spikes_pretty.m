function random_spikes_pretty(all_spikes,name,labels,montage,edf_path,out_path,mT,titles)

% seed rng
rng(0)


nspikes = 25;
surround = 4; % 2 seconds befor or after
ntotal = size(all_spikes,1);
if ntotal == 0
    fprintf('\nNo spikes detected for %s\n',name)
    return
else
    fprintf('\n%d spikes detected for %s\n',ntotal,name)
end
figure
set(gcf,'position',[10 10 1600 1000])
nwidth = 5;
nheight = 5;
t = tiledlayout(nheight,nwidth,'tilespacing','tight','padding','tight');

for i = 1:nspikes
    
    % pick a random spike
    s = randi(ntotal);

    f = all_spikes(s,4);
    ch = all_spikes(s,1);
    idx = all_spikes(s,2);
   

    % load the appropriate edf file
    path = [edf_path,name,'/file',sprintf('%d.edf',f)];
    data = edfread(path);
    info = edfinfo(path);
    samples_per_time = info.NumSamples(1,1);
    num_samples = samples_per_time*info.NumDataRecords;
    fs = round(info.NumSamples(1,1)/seconds(info.DataRecordDuration));
    rlabels = cellstr(info.SignalLabels);

    

    % get allowable electrodes
    allowable_labels = get_allowable_elecs(name);

    %% Find which of these electrodes are outside the brain
    % Remove some potentially allowable labels if they aren't really targeting mesial temporal region for that patient
    exr = strcmp(mT.name,name); assert(sum(exr==1));
    exc = mT.exclude{exr};
    if ~isempty(exc)
        C = strsplit(exc,', ');
        rm_allow = zeros(length(allowable_labels),1);
        for j = 1:length(C) % loop over exclusion labels
            rm_allow(contains(allowable_labels,C{j},'ignorecase',true)) = 1;
        end
        allowable_labels(rm_allow==1) = [];
    end

    % Turn non-A, B, C MT electrodes into A, B, C
    old_labels = rlabels;
    rlabels = mt_name_conversion(rlabels,name);

    % Find labels that match allowable electrodes and have symmetric coverage
    [allowed_labels,final_allowed_idx] = find_mt_symmetric_coverage(rlabels,allowable_labels);
    nallowed = length(allowed_labels);
    old_allowed = old_labels(final_allowed_idx);

    % Find labels that match allowable electrodes
    %allowed = ismember(rlabels,allowable_labels);
    %allowed_labels = rlabels(allowed);

    

    % Initialize values
    values = nan(num_samples,nallowed);

    % Get allowed label values
    for is = 1:nallowed
        %curr_signal = allowed_labels{is};
        curr_signal = old_allowed{is};
  
        % Fix for hup119
        if strcmp(name,'HUP119')
            if ~ismember(curr_signal,data.Properties.VariableNames)
                matching_idx = contains(data.Properties.VariableNames,curr_signal);
                curr_signal = data.Properties.VariableNames{matching_idx};
            end
        end

        %% Fill up values
        values(:,is) = data.(curr_signal){1};
        
    end

    %% Downsample to 256 hz
    if abs(fs-256)>1
        %old_values = values; old_fs = fs; old_times = linspace(0,size(old_values,1)/old_fs,size(old_values,1));
        p = 256; q = fs;
        values = resample(values,p,q);
        fs = 256;
    end

    % Convert MUSC labels into HUP labels for consistency
    allowed_labels = convert_musc_labels_to_hup(allowed_labels);

    % Get the channel index in this realigned data
    rch = find(strcmp(allowed_labels,labels{ch}));
    %rch = find(strcmp(old_allowed,labels{ch}));
    
    % apply the montage
    switch montage
        case 'machine'
            % no change
            mlabels = cellfun(@(x) sprintf('%s-Ref',x),allowed_labels,'uniformoutput',false);
        case 'car'
            [values,mlabels] = car_montage(values,1:nallowed,allowed_labels);
        case 'bipolar'
            [values,~,mlabels] = bipolar_montage_fc(values,allowed_labels,[],[],name);
    end
   
    % get the indices to plot
    start_idx = idx - fs*surround;
    end_idx = idx + fs*surround;

    times = (1:end_idx-start_idx+1)/fs;

    curr_vals = values(start_idx:end_idx,rch);
    
    %subtract median
    curr_vals = curr_vals - nanmedian(curr_vals);

    % plot the appropriate indices
    nexttile
    plot(times,curr_vals);
    hold on
    %plot(fs*surround+1,values(fs*surround+1,rch),'o')
    %title(sprintf('Spike %d %s',s,mlabels{rch}))
    xlim([0 8])
    if i > nwidth * (nheight-1)
        xlabel('Time (s)')
    else
        xticklabels([])
    end

    set(gca,'fontsize',15)

    %
    if mod(i,nwidth) == 1
        ylabel('uV')
    else
        %yticklabels([])
    end
    %}

end
title(t,titles)

print(gcf,[out_path,name,'/',montage,'_pretty'],'-dpng')
close gcf


end