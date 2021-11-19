function quick_raster(summ)

% Get the spikes and the labels
times = summ.times; times = times/3600/24; % convert times to days
spikes = summ.spikes;
labels = summ.labels;
run_dur = 60; % 60 s
spikes = spikes/run_dur*60; % convert spikes to spikes/minute (note this divides by 60 and then multiplies by 60 and so does nothing)
ekg = find_non_intracranial(labels);
sz_times = summ.sz_times/3600/24;
name = summ.name;

% remove non intracranial
labels(ekg) = [];
spikes(ekg,:) = [];
spikes = spikes/length(labels); % get # spikes/elec

set(gcf,'position',[10 10 1200 900])
    tiledlayout(1,1,'padding','tight')
    nexttile
    h = turn_nans_gray(spikes);
    hold on
    set(h,'XData',times);
    
    % show the seizure times
    for s = 1:size(sz_times,1)
        plot([sz_times(s,1) sz_times(s,1)],ylim,'r--')
        plot([sz_times(s,2) sz_times(s,2)],ylim,'r--')
    end
    
    xlim([times(1) times(end)])
    xlabel('Day')
    yticks(1:length(labels));
    yticklabels(labels);
    ylabel('Electrode')
    c = colorbar;
    ylabel(c,'Spikes/elec/min','fontsize',15);
    title(name)
    set(gca,'fontsize',15);

end