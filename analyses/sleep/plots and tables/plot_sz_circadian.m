function plot_sz_circadian

%% Parameters
plot_type = 'scatter';
nblocks = 6;
myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];



locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

%% Load out file and get roc stuff
out = load([out_folder,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];


figure

nexttile
periods = sz_circ_out.periods;
iqr_psd = sz_circ_out.iqr_psd;
all_psd = sz_circ_out.all_psd;
median_psd = nanmedian(all_psd,1);
mp = shaded_error_bars(periods,median_psd,iqr_psd,[0 0 0]);

nexttile
pre_wake = sz_circ_out.pre_wake;
n_sleep_wake = bin_out.n_sleep_wake;
perc_sz_asleep = cellfun(@(x) prc_asleep(x),pre_wake);
perc_all_asleep = 100*n_sleep_wake(:,2)./sum(n_sleep_wake,2);
minp = min([perc_sz_asleep;perc_all_asleep]);
maxp = max([perc_sz_asleep;perc_all_asleep]);
plot(perc_all_asleep,perc_sz_asleep,'o','linewidth',2)
hold on
plot([0 100],[0 100],'k--')
xlabel('Percentage of total time awake')
ylabel('Percentage of seizures arising from wakefulness')
set(gca,'fontsize',15)

nexttile
loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
p = ranksum(perc_sz_asleep(temporal),perc_sz_asleep(extra));
plot(1+randn(sum(temporal),1)*0.05,perc_sz_asleep(temporal),'o','linewidth',2,'color',myColours(1,:))
hold on
plot([0.7 1.3],[nanmedian(perc_sz_asleep(temporal)) nanmedian(perc_sz_asleep(temporal))],...
    'linewidth',2,'color',myColours(1,:))
plot(2+randn(sum(extra),1)*0.05,perc_sz_asleep(extra),'o','linewidth',2,'color',myColours(2,:))
plot([1.7 2.3],[nanmedian(perc_sz_asleep(extra)) nanmedian(perc_sz_asleep(extra))],...
    'linewidth',2,'color',myColours(2,:))
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})
ylabel('Seizures arising from wake (%)')
title('Percentage of seizures arising from wake')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)


end

function prc = prc_asleep(x)

prc = 100 * sum(x==1)/(sum(x==1)+sum(x==0));

end