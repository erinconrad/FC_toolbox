function plot_random_spikes(all_spikes,name,labels,edf_path)

nspikes = 50;
surround = 2; % 2 seconds befor or after
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
    fs = round(info.NumSamples(1,1)/seconds(info.DataRecordDuration));

    % load data from channel
    vals = data.(labels{ch});
    vals= vals{1};

    % get the indices to plot
    start_idx = idx - fs*surround;
    end_idx = idx + fs*surround;

    % plot the appropriate indices
    nexttile
    plot(vals(start_idx:end_idx));
    hold on
    plot(fs*surround+1,vals(fs*surround+1),'o')
    title(sprintf('Spike %d %s',s,labels{ch}))

end

print(gcf,[edf_path,name,'/random_spikes'],'-dpng')
close gcf


end