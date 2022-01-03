function sleep_figure4

%% Seed rng (for splitting training and testing data)
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
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
%out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

figure
set(gcf,'position',[100 100 1000 700])
tiledlayout(3,2,'tilespacing','compact','padding','compact')

%% Get stuff
rate_sw = out.bin_out.all_elecs_rates_sw;
rl_sw = out.bin_out.all_elecs_rl_sw;
rate = out.bin_out.all_elecs_rates;
rl = out.bin_out.all_elecs_rl;
leader = out.bin_out.all_elecs_leader;
%soz_rank_sw_rate = out.bin_out.soz_rank_sw;
%soz_rank_sw_rl = out.bin_out.soz_rank_sw_rl;
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

% For sleep/wake ranking comparison, generate set in which any nans in
% sleep or wake are removed (to compare same set of electrodes when
% comparing ranks in wake and asleep)
any_nans = cellfun(@(x) any(isnan(x),2),rate_sw,'uniformoutput',false);
npts = length(any_nans);
rate_sw_nan_removed = rate_sw;
rl_sw_nan_removed = rl_sw;
soz_nan_removed = soz;
rate_nan_removed = rate;
for i = 1:npts
    curr_nans = any_nans{i};
    curr_rate_sw = rate_sw{i};
    curr_rl_sw = rl_sw{i};
    curr_soz = soz{i};
    
    % remove nans
    rate_sw_nan_removed{i} = curr_rate_sw(~curr_nans,:);
    rl_sw_nan_removed{i} = curr_rl_sw(~curr_nans,:);
    soz_nan_removed{i} = curr_soz(~curr_nans);
    rate_nan_removed{i} = curr_soz(~curr_nans);
end

assert(sum(cellfun(@(x) any(isnan(x),'all'),rate_sw_nan_removed)) == 0);

wake_rate = cellfun(@(x) x(:,1), rate_sw_nan_removed,'uniformoutput',false);
sleep_rate = cellfun(@(x) x(:,2), rate_sw_nan_removed,'uniformoutput',false);


[~,wake_soz_ranks,wake_chance] = get_ranks(wake_rate,soz_nan_removed,rate_nan_removed,'rate',min_rate);
wake_soz_ranks = cellfun(@nanmedian,wake_soz_ranks);
[~,sleep_soz_ranks,sleep_chance] = get_ranks(sleep_rate,soz_nan_removed,rate_nan_removed,'rate',min_rate);
sleep_soz_ranks = cellfun(@nanmedian,sleep_soz_ranks);
soz_rank_sw_rate = [wake_soz_ranks,sleep_soz_ranks];
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
% For sleep/wake ranking comparison, generate set in which any nans in
% sleep or wake RL are also removed (to compare same set of electrodes when
% comparing ranks in wake and asleep). This removes electrodes in addition
% to the rate removal step above.
any_nans = cellfun(@(x) any(isnan(x),2),rl_sw_nan_removed,'uniformoutput',false);
npts = length(any_nans);

for i = 1:npts
    curr_nans = any_nans{i};
    curr_rate_sw = rate_sw_nan_removed{i};
    curr_rl_sw = rl_sw_nan_removed{i};
    curr_soz = soz_nan_removed{i};
    
    % remove nans
    rate_sw_nan_removed{i} = curr_rate_sw(~curr_nans,:);
    rl_sw_nan_removed{i} = curr_rl_sw(~curr_nans,:);
    soz_nan_removed{i} = curr_soz(~curr_nans);
    rate_nan_removed{i} = curr_soz(~curr_nans);
end

assert(sum(cellfun(@(x) any(isnan(x),'all'),rl_sw_nan_removed)) == 0);
assert(sum(cellfun(@(x) any(isnan(x),'all'),rate_sw_nan_removed)) == 0);

wake_rl = cellfun(@(x) x(:,1), rl_sw_nan_removed,'uniformoutput',false);
sleep_rl = cellfun(@(x) x(:,2), rl_sw_nan_removed,'uniformoutput',false);
[~,wake_soz_ranks] = get_ranks(wake_rl,soz_nan_removed,rate_nan_removed,'rl',min_rate);
wake_soz_ranks = cellfun(@nanmedian,wake_soz_ranks);
[~,sleep_soz_ranks] = get_ranks(sleep_rl,soz_nan_removed,rate_nan_removed,'rl',min_rate);
sleep_soz_ranks = cellfun(@nanmedian,sleep_soz_ranks);
soz_rank_sw_rl = [wake_soz_ranks,sleep_soz_ranks];
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
