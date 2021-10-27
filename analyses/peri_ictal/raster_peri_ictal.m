function raster_peri_ictal(out,f,s,m)

skip = 4;

data = out.file(f).sz(s).montage(m);
spikes = data.spikes;
net = data.net;
ad = data.ad;
labels = data.labels;


% Get node strength
net_uw = wrap_or_unwrap_adjacency_fc_toolbox(net);
ns = squeeze(nanmean(net_uw,1));

nruns = size(spikes,2);
sz = nruns/2;

% remove intracranial
ekg = find_non_intracranial(labels);
ns(ekg,:) = [];
spikes(ekg,:) = [];
ad(ekg,:) = [];
labels(ekg) = [];

figure
set(gcf,'position',[10 10 1200 1000])
tiledlayout(3,1,'tilespacing','compact','padding','tight')

% spikes
nexttile
turn_nans_gray(spikes)
hold on
plot([sz sz],ylim,'r--','linewidth',3)
xlabel('Minutes')
ylabel('Electrode')
title('Spike rates')
set(gca,'fontsize',15)
c1 = colorbar;
ylabel(c1,'Spikes/elec/min','fontsize',15)
yticks(1:skip:size(spikes,1))
yticklabels(labels(1:skip:size(spikes,1)))


% ns
nexttile
turn_nans_gray(ns)
hold on
plot([sz sz],ylim,'r--','linewidth',3)
xlabel('Minutes')
ylabel('Electrode')
title('Pearson correlation node strength')
set(gca,'fontsize',15)
c2 = colorbar;
ylabel(c2,'Node strength','fontsize',15)
yticks(1:skip:size(spikes,1))
yticklabels(labels(1:skip:size(spikes,1)))

% ad
nexttile
turn_nans_gray(ad)
hold on
plot([sz sz],ylim,'r--','linewidth',3)
xlabel('Minutes')
ylabel('Electrode')
title('Alpha delta ratio (lower suggests sleep or post-ictal slowing)')
set(gca,'fontsize',15)
c3 = colorbar;
ylabel(c3,'Alpha delta ratio','fontsize',15)
yticks(1:skip:size(spikes,1))
yticklabels(labels(1:skip:size(spikes,1)))


end