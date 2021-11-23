function raster_peri_ictal(pc,f,s,do_save)

%% Parameters
skip = 6; % show every skip electrode labels
m = 2; % montage (CAR)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/peri-ictal/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Put runs together
out = peri_ictal_grouping(pc);
name = pc.name;
fname = pc.file(f).name;
data = out.file(f).sz(s).montage(m);
spikes = data.spikes;
net = data.net;
ad = data.ad;
labels = data.labels;
nruns = size(spikes,2);
sz = nruns/2;
run_times = out.file(f).sz(s).run_times;

% Get all seizure times
all_szs = nan(length(pc.file(f).sz),2);
for is = 1:length(pc.file(f).sz)
    sz_start = pc.file(f).sz(is).run(nruns/2).run_times(2); % end of middle run is start of sz
    sz_end = pc.file(f).sz(is).run(nruns/2+1).run_times(1); % first of middle +1 run is end of sz
    if isnan(sz_end), sz_end = sz_start; end
    all_szs(is,:) = [sz_start sz_end];
end

% convert sz times and run centers to minutes
all_szs = all_szs/60;
run_times = run_times/60;
old_all_szs = all_szs;

% realign middle to be current sz 
all_szs = all_szs + repmat(sz-old_all_szs(s,1),size(all_szs,1),1);
% just take start
sz_to_plot = all_szs(all_szs(:,1) > 0 & all_szs(:,1) < nruns,:);

% realign later seizure based on current sz length
later_szs = sz_to_plot(:,1) > 120.1;
curr_sz = abs(sz_to_plot(:,1)-120) < 0.01;
curr_sz_length = sz_to_plot(curr_sz,2) - sz_to_plot(curr_sz,1);
%sz_to_plot(later_szs,:) = sz_to_plot(later_szs,:) + curr_sz_length;

run_times = run_times + repmat(sz-old_all_szs(s,1),size(run_times,1),2);

% confirm that one is close to center
if ~any(abs(sz_to_plot(:,1)-sz)<1), error('what'); end

% Get node strength
net_uw = wrap_or_unwrap_adjacency_fc_toolbox(net);
ns = squeeze(nanmean(net_uw,1));




% remove intracranial
ekg = find_non_intracranial(labels);
ns(ekg,:) = [];
spikes(ekg,:) = [];
ad(ekg,:) = [];
labels(ekg) = [];

figure
times = run_times(:,1);
set(gcf,'position',[10 10 1200 1000])
tiledlayout(3,1,'tilespacing','compact','padding','tight')

% spikes
nexttile
h = turn_nans_gray(spikes);
hold on
%set(h,'XData',times);
%{
for is = 1:size(sz_to_plot,1)
    plot([sz_to_plot(is,1) sz_to_plot(is,1)],ylim,'r--','linewidth',1)
    plot([sz_to_plot(is,2) sz_to_plot(is,2)],ylim,'r--','linewidth',1)
end
%}
%xlabel('Minutes')

sz_blocks = get_sz_blocks(times,sz_to_plot);
yl = ylim;
for is = 1:size(sz_to_plot,1)
    area([sz_blocks(is,1) sz_blocks(is,2)],[yl(2) yl(2)],'basevalue',yl(1),...
        'facecolor','r')
end
xlabel('Minute blocks')
ylabel(sprintf('Electrode (show every %d)',skip))
title('Spike rates')
set(gca,'fontsize',15)
c1 = colorbar;
ylabel(c1,'Spikes/elec/min','fontsize',15)
yticks(1:skip:size(spikes,1))
yticklabels(labels(1:skip:size(spikes,1)))


% ns
nexttile
h = turn_nans_gray(ns);
hold on
%set(h,'XData',times);
%{
for is = 1:size(sz_to_plot,1)
    plot([sz_to_plot(is,1) sz_to_plot(is,1)],ylim,'r--','linewidth',1)
    plot([sz_to_plot(is,2) sz_to_plot(is,2)],ylim,'r--','linewidth',1)
end
%}
%xlabel('Minutes')
sz_blocks = get_sz_blocks(times,sz_to_plot);
yl = ylim;
for is = 1:size(sz_to_plot,1)
    area([sz_blocks(is,1) sz_blocks(is,2)],[yl(2) yl(2)],'basevalue',yl(1),...
        'facecolor','r')
end

ylabel(sprintf('Electrode (show every %d)',skip))
xlabel('Minute blocks')
title('Pearson correlation node strength')
set(gca,'fontsize',15)
c2 = colorbar;
ylabel(c2,'Node strength','fontsize',15)
yticks(1:skip:size(spikes,1))
yticklabels(labels(1:skip:size(spikes,1)))

% ad
nexttile
h = turn_nans_gray(ad);
hold on
%set(h,'XData',times);
%{
for is = 1:size(sz_to_plot,1)
    plot([sz_to_plot(is,1) sz_to_plot(is,1)],ylim,'r--','linewidth',1)
    plot([sz_to_plot(is,2) sz_to_plot(is,2)],ylim,'r--','linewidth',1)
end
%}
%xlabel('Minutes')
xlabel('Minute blocks')
sz_blocks = get_sz_blocks(times,sz_to_plot);
yl = ylim;
for is = 1:size(sz_to_plot,1)
    area([sz_blocks(is,1) sz_blocks(is,2)],[yl(2) yl(2)],'basevalue',yl(1),...
        'facecolor','r')
end
ylabel(sprintf('Electrode (show every %d)',skip))
title('Alpha delta ratio (lower suggests sleep or pathological slowing)')
set(gca,'fontsize',15)
c3 = colorbar;
ylabel(c3,'Alpha delta ratio','fontsize',15)
yticks(1:skip:size(spikes,1))
yticklabels(labels(1:skip:size(spikes,1)))

if do_save
    print(gcf,[out_folder,sprintf('%s_sz%d',fname,s)],'-dpng');
    close(gcf)
end

end