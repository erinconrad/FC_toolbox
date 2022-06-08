function do_model(from_scratch,just_gray)
%{



%}

rng(0) % seed rng (for bootstrap portion)

%% Parameters

durations = {1, 5, 10, 30, 60,[]};
ndurs = length(durations);
ncoeffs = 4;
myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];



%% Locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];
time_roc_folder = [out_folder,'time_roc/'];
if ~exist(time_roc_folder,'dir')
    mkdir(time_roc_folder)
end

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;
npts = length(out.circ_out.names);
all_soz = out.bin_out.all_is_soz;

%% Get number of patients with electrode localizations
has_elec_locs = nan(npts,1);
for ip = 1:npts
    curr_locs = out.circ_out.all_elec_locs{ip};
    not_all_empty = any(~cellfun(@isempty,curr_locs));
    has_elec_locs(ip) = not_all_empty;
end


%% Get some clinical stuff
stereo = logical(out.circ_out.stereo);
loc = out.circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');

%% Do the stats
if from_scratch
    % Bootstrap stats
    fprintf('\nDoing bootstrap to get model stats\n');
    [coeff_names,coeff_stats] = sleep_model_bootstrap_stats(just_gray);

    % LOO
    fprintf('\nDoing LOO analysis to get individual patient performance\n');
    [pt_stats,X,Y,pt_specific,excluded] = sleep_loo(just_gray);

    % Get AUC for sleep and wake as a function of duration
    fprintf('\nDoing duration analysis\n');
    time_aucs = sleep_duration(durations,just_gray);

    nmout.coeff_stats = coeff_stats;
    nmout.coeff_names = coeff_names;
    nmout.pt_stats = pt_stats;
    nmout.X = X;
    nmout.Y = Y;
    nmout.pt_specific = pt_specific;
    nmout.time_aucs = time_aucs;
    nmout.excluded = excluded;

    if just_gray
        save([time_roc_folder,'new_model_out_gray.mat'],'nmout');
    else
        save([time_roc_folder,'new_model_out.mat'],'nmout');
    end
    
else
    
    nmout_gray = load([time_roc_folder,'new_model_out_gray.mat']);
    nmout = load([time_roc_folder,'new_model_out.mat']);
    nmout = nmout.nmout;
    nmout_gray = nmout_gray.nmout;
    unpack_any_struct(nmout);
end

if from_scratch
    return
end

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');
sfid = fopen([out_folder,'supplemental_results.html'],'a');
fprintf(fid,'<p><b>Localizing the seizure onset zone with spikes</b><br>');

%% Print results of main model
fprintf(fid,['We tested how accurately spike rates could classify electrodes '...
    'as SOZ versus non-SOZ. To do this, we constructed a mixed-effects logistic regression classifier, '...
    'where each electrode from each patient was a separate observation. The response '...
    'was a binary variable indicating whether the electrode was SOZ or non-SOZ,'...
    ' and the predictors were the average spike rate in wakefulness, the average spike '...
    'rate in sleep, the average preictal (6 hour window) spike rate, and the average '...
    'postictal spike rate. Spike rates were normalized across electrodes within each patient. '...
    'The patient identifier was used as a random effect. To calculate statistics on model coefficients, we used bootstrapping (N = 1,000), '...
    'where for each bootstrap sample we randomly selected patients with replacement. '...
    'Holding other predictors constant, the odds of an electrode being a SOZ increased significantly '...
    'with an increase in spike rate in sleep (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s) '...
    'and with an increase in postictal spike rate (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s), '...
    'but not with an increase in spike rate in wakefulness (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s) '...
    'or an increase in preictal spike rate (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s). These results imply '...
    'that spikes in sleep and postictal spikes independently localize the SOZ, but spikes in wakefulness '...
    'and preictal spikes do not. When restricting analysis to only electrodes in gray matter, '...
    'we again found that spike rates in sleep were associated with a higher likelihood of being '...
    'in the SOZ, but we found no association for spikes in any other state (Supplemental Results).</p>'], ...
    (exp(coeff_stats(1,1))-1)*100,exp(coeff_stats(1,2)),exp(coeff_stats(1,3)),get_p_html(coeff_stats(1,4)),...
    (exp(coeff_stats(4,1))-1)*100,exp(coeff_stats(4,2)),exp(coeff_stats(4,3)),get_p_html(coeff_stats(4,4)),...
    (exp(coeff_stats(2,1))-1)*100,exp(coeff_stats(2,2)),exp(coeff_stats(2,3)),get_p_html(coeff_stats(2,4)),...
    (exp(coeff_stats(3,1))-1)*100,exp(coeff_stats(3,2)),exp(coeff_stats(3,3)),get_p_html(coeff_stats(3,4)));


%% Print results of gray matter analysis
fprintf(sfid,'<p><b>Localizing the seizure onset zone with spikes - gray matter only</b><br>');
fprintf(sfid,['We repeated the SOZ localization analysis, only including electrode '...
    'contacts determined to be in gray matter from clinical imaging reconstruction. '...
    'Holding other predictors constant, the odds of an electrode being a SOZ increased significantly '...
    'with an increase in spike rate in sleep (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s) '...
    'and with an increase in postictal spike rate (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s), '...
    'but not with an increase in spike rate in wakefulness (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s) '...
    'or an increase in preictal spike rate (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s).</p>'], ...
    (exp(nmout_gray.coeff_stats(1,1))-1)*100,exp(nmout_gray.coeff_stats(1,2)),exp(nmout_gray.coeff_stats(1,3)),get_p_html(nmout_gray.coeff_stats(1,4)),...
    (exp(nmout_gray.coeff_stats(4,1))-1)*100,exp(nmout_gray.coeff_stats(4,2)),exp(nmout_gray.coeff_stats(4,3)),get_p_html(nmout_gray.coeff_stats(4,4)),...
    (exp(nmout_gray.coeff_stats(2,1))-1)*100,exp(nmout_gray.coeff_stats(2,2)),exp(nmout_gray.coeff_stats(2,3)),get_p_html(nmout_gray.coeff_stats(2,4)),...
    (exp(nmout_gray.coeff_stats(3,1))-1)*100,exp(nmout_gray.coeff_stats(3,2)),exp(nmout_gray.coeff_stats(3,3)),get_p_html(nmout_gray.coeff_stats(3,4)));

%% initialize figure
figure
set(gcf,'position',[10 10 1000 1200])
tiledlayout(3,2,'tilespacing','tight','padding','tight')


%% LOO
aucs = pt_stats(:,5);
scores =  pt_specific(:,1);
soz = pt_specific(:,3);
threshold = pt_specific(:,2);
desired_threshold = [];
[npv,ppv,pred,totaln,mat,threshold,acc] = cellfun(@(x,y,z) individual_threshold_stats(x,y,z,desired_threshold),...
    scores,soz,threshold,'uniformoutput',false);
npv = cell2mat(npv);
ppv = cell2mat(ppv);
acc = cell2mat(acc);
threshold_mat = cell2mat(threshold);

%% PLot ROC
ym = nanmean(Y,1);
ly = ym-nanstd(Y,[],1);
uy = ym+nanstd(Y,[],1);
nexttile
%plot(X,ym,'linewidth',2)
mp = shaded_error_bars_fc(X,ym,[ly;uy],'k');
hold on
plot([0 1],[0 1],'k--','linewidth',2)
legend(mp,sprintf('Mean AUC %1.2f',nanmean(aucs)),'location','southeast','fontsize',15)
xlabel('False positive rate')
ylabel('True positive rate')
title('SOZ classification accuracy')
set(gca,'fontsize',15)

%% Text
fprintf(fid,['<p>We next assessed how a model well a incorporating spike rates would '...
    'localize the SOZ in new patients. We used a leave-one-out training model where, for '...
    'each patient, we trained the logistic regression classifier described above on all other patients '...
    'and tested the classifier on the held-out patient. We excluded patients without seizures '...
    '(N = %d), patients without detected sleep (N = %d), and patients for whom the model '...
    'failed to converge (N = %d), resulting in %d remaining patients. The mean (standard deviation) AUC of '...
    'an ROC curve representing classifier performance on the test data was %1.2f (%1.2f) (Fig 5A). '...
    'Results were similar for a model analyzing electrodes in gray matter only '...
    '(mean (SD) AUC =  %1.2f (%1.2f)).</p>'],...
    4,1,1,sum(~isnan(aucs)),nanmean(aucs,1),nanstd(aucs,[],1),nanmean(nmout_gray.pt_stats(:,5),1),nanstd(nmout_gray.pt_stats(:,5),[],1));

%% Get all data for single confusion matrix
% Convert everything to long vectors
all_scores = [];
all_thresholds = [];
all_soz = [];
for i = 1:length(soz)
    all_scores = [all_scores;scores{i}];
    all_soz = [all_soz;soz{i}];
    all_thresholds = [all_thresholds;repmat(threshold{i},length(soz{i}),1)];
    
end

% Do confusion matrix for this
[all_npv,all_ppv,all_pred,all_totaln,all_mat,all_threshold,all_acc] = individual_threshold_stats(all_scores,all_soz,[],all_thresholds);

% Compare to alternate way, make sure same
mat = cat(3,mat{:});
mat = nansum(mat,3);
alt_total_ppv = mat(2,2)/(mat(2,2) + mat(1,2));
alt_total_npv = mat(1,1)/(mat(1,1) + mat(2,1));

assert(alt_total_ppv == all_ppv || alt_total_npv == all_npv)

%% Plot confusion matrix
nexttile
turn_nans_gray([1 0;0 1])
colormap(gca,[0.8000, 0.420, 0.42;0.5, 0.75, 0.5])
xticks(1:2)
xticklabels({'Not SOZ','SOZ'})
yticks(1:2)
yticklabels({'Not SOZ','SOZ'})
xlabel('Predicted')
ylabel('Actual')
hold on
for ic = 1:2
    for jc = 1:2
        text(ic,jc,sprintf('%d',all_mat(jc,ic)),'horizontalalignment','center','fontsize',25,'fontweight','bold')
    end
end
title(sprintf('All patients in aggregate\nAccuracy: %1.1f%%, PPV: %1.1f%%, NPV: %1.1f%%',...
all_acc*100,...
all_ppv*100,all_npv*100))
set(gca,'fontsize',15)


%% Text
fprintf(fid,['<p>To provide an example of real-world classifier performance, we set a'...
    ' probability threshold and calculated a confusion matrix. We selected patient-specific thresholds, '...
    'each chosen such that the number of SOZ electrodes predicted by the model equaled the number of '...
    'clinician-defined SOZ electrodes (which varied across patients). This led to a mean (SD) '...
    'threshold probability of %1.2f (%1.2f) across patients. Aggregating all patient data, the '...
    'resulting confusion matrix had an accuracy of %1.1f%%, a positive predictive value (PPV) of '...
    '%1.1f%%, and a negative predictive value (NPV) of %1.1f%% (Fig 5B). Note that the large discrepancy between '...
    'PPV and NPV reflects the fact that most electrodes do not belong to the SOZ.</p>'],...
    nanmean(threshold_mat), nanstd(threshold_mat),all_acc*100,all_ppv*100,all_npv*100);



%% PPVs and NPVs for each patient
% Sort lowest to highest AUC
nexttile
%[~,I] = sort(ppv);
%plot(1:length(ppv),ppv(I),'o','linewidth',2)
%hold on
%plot(1:length(ppv),npv(I),'o','linewidth',2)
plot(1:length(aucs),sort(aucs),'ko','linewidth',2)
%legend({'PPV','NPV'},'fontsize',15,'location','southeast')
title('Individual patient accuracies')
set(gca,'fontsize',15)
xlabel('Patient')
ylabel('AUC')

fprintf(fid,'<p>Model accuracies were highly variable across patients (Fig 5C.');

%% Do between-patient analyses

% Implant type
nexttile
plot(1+randn(sum(stereo),1)*0.05,aucs(stereo),'o','linewidth',2,'color',myColours(1,:))
hold on
plot(2+randn(sum(~stereo),1)*0.05,aucs(~stereo),'o','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
xticks([1 2])
xticklabels({'Stereo-EEG','Grid/strip/depths'});
ylabel('Classification AUC')
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
nylim = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k','linewidth',2)
text(1.5,ytext,get_p_text(ranksum(aucs(stereo),aucs(~stereo))),...
    'fontsize',15,'horizontalalignment','center')
ylim(nylim)
xlim([0.5 2.5])
title('Accuracy by implant strategy')
set(gca,'fontsize',15)

[p,~,stats] = ranksum(aucs(stereo),aucs(~stereo));
% text
W = stats.ranksum;
nt = sum(~(isnan(aucs(stereo))));
ne = sum(~(isnan(aucs(~stereo))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,[' Patients with stereo-EEG implantations had similar '...
    'classification accuracy '...
    '(median AUC = %1.2f) to those with grid/strip/depth implantations'...
    ' (median AUC = %1.2f) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>stereo</sub></i> = %d, <i>N<sub>not-stereo</sub></i> = %d) ='...
    ' %1.1f, %s) (Fig 5D).'],nanmedian(aucs(stereo)),nanmedian(aucs(~stereo)),...
    nt,ne,U,get_p_html(p));

% Localization
nexttile
plot(1+randn(sum(temporal),1)*0.05,aucs(temporal),'o','linewidth',2,'color',myColours(2,:))
hold on
plot(2+randn(sum(extra),1)*0.05,aucs(extra),'o','linewidth',2,'color',myColours(3,:))
xticks([1 2])
xticklabels({'Temporal lobe epilepsy','Extra-temporal lobe epilepsy'});
ylabel('Classification AUC')
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
nylim = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k','linewidth',2)
text(1.5,ytext,get_p_text(ranksum(aucs(temporal),aucs(extra))),...
    'fontsize',15,'horizontalalignment','center')
ylim(nylim)
xlim([0.5 2.5])
title('Accuracy by seizure localization')
set(gca,'fontsize',15)
[p,~,stats] = ranksum(aucs(temporal),aucs(extra));
% text
W = stats.ranksum;
nt = sum(~(isnan(aucs(temporal))));
ne = sum(~(isnan(aucs(extra))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,[' Patients with temporal lobe epilepsy had higher '...
    'classification accuracy '...
    '(median AUC = %1.2f) relative to those with extra-temporal'...
    ' lobe epilepsy (median AUC = %1.2f) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>temporal</sub></i> = %d, <i>N<sub>extra-temporal</sub></i> = %d) ='...
    ' %1.1f, %s) (Fig 5E).</p>'],nanmedian(aucs(temporal)),nanmedian(aucs(extra)),...
    nt,ne,U,get_p_html(p));


%% Duration
% Prep duration ticks
duration_ticks = cell2mat(durations(1:ndurs-1));
new_lim = 1.2*(max(duration_ticks)-0)+min(duration_ticks);
duration_ticks = [duration_ticks,new_lim];
duration_ticklabels = arrayfun(@(x) sprintf('%d',x),duration_ticks,'uniformoutput',false);
duration_ticklabels(end) = {'Full'};

nexttile
sp = nan(2,1);
for iws = 1:2
    curr_aucs = squeeze(nanmean(time_aucs(iws,:,:,:),[3 4]));
    auc_std = squeeze(nanstd(time_aucs(iws,:,:,:),[],[3 4]));
    sp(iws) = plot(duration_ticks,curr_aucs,'-o','linewidth',2);
    hold on
end
legend({'Wake','Sleep'},'fontsize',15,'location','southeast')
xticks(duration_ticks)
xticklabels(duration_ticklabels)
ylabel('Mean AUC')
xlabel('Duration examined (minutes)')
title('Accuracy by duration')
set(gca,'fontsize',15);

% Text
%
fprintf(fid,['<p>We next studied what duration of interictal data was needed '...
    'to localize the SOZ. Given that spikes in sleep better localized the SOZ than spikes in wakefulness, '...
    'we separately evaluated sleep and wake interictal periods. We performed the same logistic regression '...
    'classification described above, again using a leave-one-out approach to train and test the model.'...
    ' Using durations of 1 minute, 5 minutes, 10 minutes, and 60 minutes, we chose random 1-minute segments '...
    'of interictal data of the specified state (sleep or wake) adding up to the specified duration. We performed '...
    'each training/testing split ten times to obtain a distribution of random durations. '...
    'We measured the AUC of the model tested on the testing data, and averaged this against each of the '...
    '10 random samples and each of the training/testing splits. This resulted in a single '...
    'average model AUC per duration and per wake vs. sleep state. The average AUC '...
    'was higher at each duration for models trained on sleep spike data'...
    ' than those trained on wake spike data. Also, near-optimal model performance was '...
    'achieved with only 10 minutes of interictal data (mean (SD) AUC for model trained on '...
    'full duration of sleep data = %1.2f (%1.2f) versus ten minutes of sleep '...
    'data = %1.2f (%1.2f)) (Fig 5F).</p>'],...
    nanmean(time_aucs(2,end,:,:),[3 4]),nanstd(time_aucs(2,end,:,:),[],[3 4]),...
    nanmean(time_aucs(2,3,:,:),[3 4]),nanstd(time_aucs(2,3,:,:),[],[3 4]));
%}

%% Annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.9 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.55 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.55 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.23 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.23 0.1 0.1],'String','F','fontsize',25,'linestyle','none')


%% Print figure
if just_gray
    print(gcf,[time_roc_folder,'Fig5'],'-dpng')

else
    print(gcf,[time_roc_folder,'Fig5'],'-dpng')
end

close all

fclose(fid);
fclose(sfid);

end