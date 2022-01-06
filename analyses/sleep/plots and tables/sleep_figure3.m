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

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');
fprintf(fid,'<p><b>Localizing the seizure onset zone with spikes</b><br>');

figure
set(gcf,'position',[100 100 800 700])
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
nexttile([1 1])
stats = plot_paired_data(soz_rate_sw',{'wake','sleep'},sprintf('Spikes/elecs/min'),'paired','scatter');
title(sprintf('SOZ spike rate'))

% Results text
fprintf(fid,['We next compared the ability of spikes to localize the SOZ '...
    'in wake versus in sleep. The spike rate of SOZ electrodes was higher '...
    ' in sleep (median %1.1f spikes/elecs/min)'...
    ' than wake (median %1.1f spikes/elecs/min) '...
    '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 4A).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

%% Plot non SOZ spike rate wake and sleep
nexttile([1 1])
stats = plot_paired_data(non_soz_rate_sw',{'wake','sleep'},sprintf('Spikes/elecs/min'),'paired','scatter');
title(sprintf('Non-SOZ spike rate'))

% Results text
fprintf(fid,[' The spike rate of non-SOZ electrodes was also higher '...
    ' in sleep (median %1.1f spikes/elecs/min)'...
    ' than wake (median %1.1f spikes/elecs/min) '...
    '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 4B).</p>'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));


%% SOZ spike rate ranking
nexttile([1 2])
stats_out = plot_orders(rate,soz);
hold on
xticklabels([])
xlabel('Patient')
ylabel('Electrode spike rate rank')
set(gca,'fontsize',15)
title('Seizure onset zone - spike rate ranking')

fprintf(fid,['<p>To determine whether sleep disproportionately increased spikes '...
    'in the SOZ relative to non-SOZ electrodes, we examined the ranking of electrodes'...
    ' by spike rates. We first examined whether SOZ electrodes '...
    'had a higher spike rate than expected by chance. We compared the median '...
    'spike rate rank of the SOZ electrodes against the median '...
    'rank of all electrodes. The median SOZ electrode '...
    'rank across patients was %1.1f, implying that across patients, %d '...
    'non-SOZ electrodes tended to have higher spike rates than the SOZ electrodes. '...
    'The median SOZ rank was higher (closer to 1) than the median overall electrode rank '...
    'in %d of %d patients, which is more than expected by chance '...
    '(Binomial test, %s) (Fig. 4C). These findings imply that the SOZ has '...
    'more frequent spikes than expected by chance.'],stats_out.median_rank,stats_out.median_rank-1,...
    stats_out.nsuc,...
    stats_out.n,get_p_html(stats_out.pval));



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
nexttile([1 1])
stats = plot_paired_data(soz_rank_sw_rate',{'wake','sleep'},sprintf('Spike rate rank'),'paired','scatter','ranking');
title(sprintf('SOZ spike rate ranking'))

% Results text
fprintf(fid,[' We next compared the ranking of SOZ electrodes '...
    'by spike rate in wake versus sleep. The spike rate ranking of SOZ electrodes was higher '...
    ' (closer to 1) in sleep (median rank %1.1f)'...
    ' than in wake (median rank %1.1f) '...
    '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 4D).'...
    ' This implies that sleep preferentially increases spike rates '...
    'in the SOZ relative to other electrodes.</p>'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

% confirm that the total number of electrodes being ranked is the same for
% wake and sleep (otherwise unfair comparison)
assert(isequal(cellfun(@length,wake_all_ranks),cellfun(@length,sleep_all_ranks)))

%% Fancy model
soz_roc_out = classify_soz;
roc = soz_roc_out.roc;
auc = soz_roc_out.auc;
sleep_or = soz_roc_out.sleep_or;
wake_or = soz_roc_out.wake_or;
sleep_ci95 = soz_roc_out.sleep_ci95;
wake_ci95 = soz_roc_out.wake_ci95;
sleep_t= soz_roc_out.sleep_t;
wake_t= soz_roc_out.wake_t;
sleep_p= soz_roc_out.sleep_p;
wake_p= soz_roc_out.wake_p;

sleep_or = (sleep_or - 1)*100;
wake_or = (wake_or - 1)*100;

nexttile([1 1])
plot(roc(:,1),roc(:,2),'k-','linewidth',2)
hold on
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','southeast','fontsize',15)
set(gca,'fontsize',15)
title('SOZ identification accuracy')

%% ROC text
fprintf(fid,['<p>Finally, we tested how accurately spike rates could '...
    'classify electrodes as SOZ versus non-SOZ. We trained a logistic regression'...
    ' classifier on a randomly chosen two-thirds of the patients, using the normalized '...
    'electrode-specific average spike rate in each of the sleep and wake states as predictor variables '...
    'and the binary classification of the electrode as SOZ or non-SOZ as '...
    'the response variable. Holding spike rate in wake constant, the odds'...
    ' of an electrode being a SOZ electrode increased by %1.1f%% (95%% CI [%1.2f-'...
    '%1.2f]) for each additional normalized spike rate in sleep unit, which was significant '...
    '(<i>t</i> = %1.1f, %s). Holding spike rate in sleep constant,'...
    ' the odds of an electrode being a SOZ electrode increased by %1.1f%% (95%% CI [%1.2f-'...
    '%1.2f]) for each additional normalized spike rate in wake unit, which was not significant ('...
    '<i>t</i> = %1.1f, %s). Applying this model to the remaining'...
    ' one-third of patients held as testing data accurately '...
    'classified electrodes as SOZ versus non-SOZ '...
    '(AUC = %1.2f, Fig. 4E). This result implies that using '...
    'spike rates and sleep/wake classification alone can accurately identify the SOZ.'...
    ' It also provides additional support that spikes in sleep particularly localize the SOZ.</p>'],...
    sleep_or,sleep_ci95(1),sleep_ci95(2),sleep_t,get_p_html(sleep_p),wake_or,...
    wake_ci95(1),wake_ci95(2),wake_t,get_p_html(wake_p),auc);

%% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0.49 0.91 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.59 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.245 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
annotation('textbox',[0.49 0.245 0.1 0.1],'String','E','fontsize',20,'linestyle','none')
%annotation('textbox',[0.49 0.24 0.1 0.1],'String','F','fontsize',20,'linestyle','none')

fclose(fid);
print([out_folder,'fig4'],'-depsc')

end
