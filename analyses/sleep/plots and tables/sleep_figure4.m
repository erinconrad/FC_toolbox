function sleep_figure4

%% Seed rng (because it's a MC test)
rng(0)
nb = 1e4;
min_rate = 0.1;


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
set(gcf,'position',[100 100 1000 700])
tiledlayout(3,2,'tilespacing','compact','padding','compact')

%% Get stuff
rate = out.bin_out.all_elecs_rates;
rl = out.bin_out.all_elecs_rl;
leader = out.bin_out.all_elecs_leader;
soz_rank_sw_rate = out.bin_out.soz_rank_sw;
soz_rank_sw_rl = out.bin_out.soz_rank_sw_rl;
soz = out.bin_out.all_is_soz;
locs = out.circ_out.all_locs;
is_temporal = cellfun(@(x) strcmp(x,'temporal'),locs);

%% SOZ spike rate ranking
nexttile([1 1])
plot_orders(rate,soz,rate,'rate',[])
hold on
xticklabels([])
xlabel('Patient')
ylabel('Electrode spike rate rank')
set(gca,'fontsize',15)
title('Seizure onset zone - spike rate ranking')

%% Spike rate ranking sleep vs wake
nexttile
plot_paired_data(soz_rank_sw_rate',{'wake','sleep'},sprintf('Spike rate rank'),'paired','scatter','ranking')
title(sprintf('SOZ spike rate ranking by wake vs sleep'))

%% SOZ spike timing ranking
nexttile([1 1])
plot_orders(rl,soz,rate,'rl',min_rate)
xticklabels([])
xlabel('Patient')
ylabel('Electrode spike timing rank')
set(gca,'fontsize',15)
title('Seizure onset zone - spike timing ranking')

%% Spike timing ranking sleep vs wake
nexttile
plot_paired_data(soz_rank_sw_rl',{'wake','sleep'},sprintf('Spike timing rank'),'paired','scatter','ranking')
title(sprintf('SOZ spike timing ranking by wake vs sleep'))

%% Correlate rate and rl
rl_rate_corr = nan(length(rate),1);
for ip = 1:length(rate)
    curr_rate = rate{ip};
    curr_rl = rl{ip};
    % Only do it for those with enough spikes
    spikey = curr_rate > min_rate;
    rl_rate_corr(ip) = corr(curr_rate(spikey),curr_rl(spikey),'type','spearman','rows','pairwise');
end

% Plot it
nexttile([1 1])
plot(rl_rate_corr,'o','linewidth',2)
hold on
plot(xlim,[0 0],'k--','linewidth',2)
ylim([-1 1])
xticklabels([])
xlabel('Patient')
ylabel('Correlation coefficient');
set(gca,'fontsize',15)
title({'Correlation between spike rate and spike timing'})



%% Fancy model
soz_roc_out = classify_soz;
roc = soz_roc_out.roc;
auc = soz_roc_out.auc;

nexttile([1 1])
plot(roc(:,1),roc(:,2),'k-','linewidth',2)
hold on
plot([0 1],[0 1],'k--','linewidth',2)
%plot(roc(disc_I,1),roc(disc_I,2),'*','markersize',15,'linewidth',2,'color',colors(5,:));
%text(roc(disc_I,1)+0.01,roc(disc_I,2)-0.05,'SOZ cutoff','fontsize',15,'color',colors(5,:));
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','southeast','fontsize',15)
set(gca,'fontsize',15)
title('SOZ identification accuracy')

%% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0.49 0.91 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.59 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
annotation('textbox',[0.49 0.59 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.24 0.1 0.1],'String','E','fontsize',20,'linestyle','none')
annotation('textbox',[0.49 0.24 0.1 0.1],'String','F','fontsize',20,'linestyle','none')

print([out_folder,'fig4'],'-dpng')

end
