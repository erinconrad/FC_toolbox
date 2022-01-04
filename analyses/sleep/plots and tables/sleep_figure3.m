function sleep_figure3



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
rate = out.bin_out.all_elecs_rates;
soz = out.bin_out.all_is_soz;
npts = length(soz);

%% Get SOZ and non SOZ spike rate wake and sleep
% initialize
soz_rate_sw = nan(npts,2);
non_soz_rate_sw = nan(npts,2);

% Loop over patients
for i = 1:npts
    curr_rate_sw = rate_sw{i};
    curr_soz = soz{i};
    soz_rate_sw(i,:) = nanmean(curr_rate_sw(curr_soz,:),1);
    non_soz_rate_sw(i,:) = nanmean(curr_rate_sw(~curr_soz,:),1);
end

%% Plot SOZ spike rate wake and sleep
nexttile
plot_paired_data(soz_rate_sw',{'wake','sleep'},sprintf('Spikes/elecs/min'),'paired','scatter')
title(sprintf('SOZ spike rate'))


%% Plot non SOZ spike rate wake and sleep
nexttile
plot_paired_data(non_soz_rate_sw',{'wake','sleep'},sprintf('Spikes/elecs/min'),'paired','scatter')
title(sprintf('Non-SOZ spike rate'))

%% SOZ spike rate ranking
nexttile([1 1])
plot_orders(rate,soz)
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
soz_nan_removed = soz;
rate_nan_removed = rate;
for i = 1:npts
    curr_nans = any_nans{i};
    curr_rate_sw = rate_sw{i};
    curr_soz = soz{i};
    
    % remove nans
    rate_sw_nan_removed{i} = curr_rate_sw(~curr_nans,:);
    soz_nan_removed{i} = curr_soz(~curr_nans);
    rate_nan_removed{i} = curr_soz(~curr_nans);
end

% confirm no nans left for wake or sleep
assert(sum(cellfun(@(x) any(isnan(x),'all'),rate_sw_nan_removed)) == 0);

% Get wake and sleep rate
wake_rate = cellfun(@(x) x(:,1), rate_sw_nan_removed,'uniformoutput',false);
sleep_rate = cellfun(@(x) x(:,2), rate_sw_nan_removed,'uniformoutput',false);


[wake_all_ranks,wake_soz_ranks] = simple_rate_rank(wake_rate,soz_nan_removed);
wake_soz_ranks = cellfun(@nanmedian,wake_soz_ranks);

[sleep_all_ranks,sleep_soz_ranks] = simple_rate_rank(sleep_rate,soz_nan_removed);
sleep_soz_ranks = cellfun(@nanmedian,sleep_soz_ranks);

soz_rank_sw_rate = [wake_soz_ranks,sleep_soz_ranks];

% Plot it
nexttile
plot_paired_data(soz_rank_sw_rate',{'wake','sleep'},sprintf('Spike rate rank'),'paired','scatter','ranking')
title(sprintf('SOZ spike rate ranking by wake vs sleep'))

% confirm that the total number of electrodes being ranked is the same for
% wake and sleep (otherwise unfair comparison)
assert(isequal(cellfun(@length,wake_all_ranks),cellfun(@length,sleep_all_ranks)))

%% Fancy model
soz_roc_out = classify_soz;
roc = soz_roc_out.roc;
auc = soz_roc_out.auc;

nexttile([1 2])
plot(roc(:,1),roc(:,2),'k-','linewidth',2)
hold on
plot([0 1],[0 1],'k--','linewidth',2)
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
%annotation('textbox',[0.49 0.24 0.1 0.1],'String','F','fontsize',20,'linestyle','none')

print([out_folder,'fig4'],'-dpng')

end
