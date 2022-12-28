function plot_random_spikes(all_spikes,name,labels,montage,edf_path,edf_summ_path)

nspikes = 50;
surround = 3; % 2 seconds befor or after
ntotal = size(all_spikes,1);
if ntotal == 0
    fprintf('\nNo spikes detected for %s\n',name)
    return
else
    fprintf('\n%d spikes detected for %s\n',ntotal,name)
end
figure
set(gcf,'position',[10 10 1400 1000])
tiledlayout(5,10,'tilespacing','tight','padding','tight')

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
    allowable_labels = get_allowable_elecs;

    % Find labels that match allowable electrodes
    allowed = ismember(rlabels,allowable_labels);
    allowed_labels = rlabels(allowed);
    nallowed = sum(allowed);

    % Initialize values
    values = nan(num_samples,nallowed);

    % Get allowed label values
    for is = 1:nallowed
        curr_signal = allowed_labels{is};
  
        %% Fill up values
        values(:,is) = data.(curr_signal){1};
        
    end

    % Get the channel index in this realigned data
    rch = find(strcmp(allowed_labels,labels{ch}));
    
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

    % plot the appropriate indices
    nexttile
    plot(values(start_idx:end_idx,rch));
    hold on
    plot(fs*surround+1,values(fs*surround+1,rch),'o')
    title(sprintf('Spike %d %s',s,mlabels{rch}))

end

print(gcf,[edf_summ_path,name,'/',montage],'-dpng')
close gcf


end