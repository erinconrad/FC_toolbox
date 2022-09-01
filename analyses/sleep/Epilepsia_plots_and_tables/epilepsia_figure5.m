function epilepsia_figure5


%% Parameters
myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];


%% Locations
locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/epilepsia/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;
npts = length(out.circ_out.names);

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

%% Get the data
nmout = out.model_out;
nmout_gray = out.model_out_gray;
    
durations = nmout.durations;
ndurs = length(durations);

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');
sfid = fopen([out_folder,'supplemental_results.html'],'a');
fprintf(fid,'<p><b>Localizing the seizure onset zone with spikes</b><br>');

%% Print results of main model
coeff_stats = nmout.coeff_stats;
fprintf(fid,['We tested how accurately spike rates could classify electrodes '...
    'as SOZ versus non-SOZ. Holding other '...
    'predictors constant, the odds of an electrode being in the SOZ increased significantly '...
    'with an increase in normalized spike rate in sleep (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s) '...
    'and with an increase in postictal spike rate (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s), '...
    'but not with an increase in spike rate in wakefulness (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s) '...
    'or an increase in preictal spike rate (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s). These results imply '...
    'that spikes in sleep and postictal spikes independently localize the SOZ, but spikes in wakefulness '...
    'and preictal spikes do not. When restricting analysis to only electrodes in gray matter, '...
    'we again found that spike rates in sleep were associated with a higher likelihood of being '...
    'in the SOZ. However, in this case spikes in the postictal '...
    'state were not associated with the likelihood of an electrode belonging to the SOZ (Supplemental Results).'], ...
    (exp(coeff_stats(1,1))-1)*100,exp(coeff_stats(1,2)),exp(coeff_stats(1,3)),get_p_html(coeff_stats(1,4)),...
    (exp(coeff_stats(4,1))-1)*100,exp(coeff_stats(4,2)),exp(coeff_stats(4,3)),get_p_html(coeff_stats(4,4)),...
    (exp(coeff_stats(2,1))-1)*100,exp(coeff_stats(2,2)),exp(coeff_stats(2,3)),get_p_html(coeff_stats(2,4)),...
    (exp(coeff_stats(3,1))-1)*100,exp(coeff_stats(3,2)),exp(coeff_stats(3,3)),get_p_html(coeff_stats(3,4)));


fprintf(fid,[' Concordant localization of preimplant data was associated with a decreased likelihood of an individual electrode '...
    'belonging to the SOZ (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s), concordant lateralization '...
    'was associated with an increased likelihood of an individual electrode belonging to the SOZ '...
    '(%1.1f%%, 95%% CI [%1.2f-%1.2f], %s), and MRI lesional status had no association (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s).'],...
    (exp(coeff_stats(6,1))-1)*100,exp(coeff_stats(6,2)),exp(coeff_stats(6,3)),get_p_html(coeff_stats(6,4)),...
    (exp(coeff_stats(7,1))-1)*100,exp(coeff_stats(7,2)),exp(coeff_stats(7,3)),get_p_html(coeff_stats(7,4)),...
    (exp(coeff_stats(5,1))-1)*100,exp(coeff_stats(5,2)),exp(coeff_stats(5,3)),get_p_html(coeff_stats(5,4)));

%% Print results of gray matter analysis
perc_in_gray =count_gray_white(out);
fprintf(sfid,'<p><b>Localizing the seizure onset zone with spikes - gray matter only</b><br>');
fprintf(sfid,['We repeated the SOZ localization analysis, only including electrode '...
    'contacts determined to be in gray matter from clinical imaging reconstruction. '...
    ' Of the 95 patients analyzed in the original analysis, 54 patients had clinical reconstructions '...
    'available indicating gray/white matter designations and were included in this secondary analysis. On average %1.1f%% of electrode contacts were '...
    'deemed to be in gray matter. '...
    'Holding other predictors constant, the odds of an electrode being a SOZ increased significantly '...
    'with an increase in spike rate in sleep (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s), '...
    'but not with an increase in spike rate in wakefulness (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s), '...
    'an increase in preictal spike rate (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s), '...
    'or an increase in postictal spike rate (%1.1f%%, 95%% CI [%1.2f-%1.2f], %s).</p>'], ...
    perc_in_gray*100,...
    (exp(nmout_gray.coeff_stats(1,1))-1)*100,exp(nmout_gray.coeff_stats(1,2)),exp(nmout_gray.coeff_stats(1,3)),get_p_html(nmout_gray.coeff_stats(1,4)),...
    (exp(nmout_gray.coeff_stats(2,1))-1)*100,exp(nmout_gray.coeff_stats(2,2)),exp(nmout_gray.coeff_stats(2,3)),get_p_html(nmout_gray.coeff_stats(2,4)),...
    (exp(nmout_gray.coeff_stats(3,1))-1)*100,exp(nmout_gray.coeff_stats(3,2)),exp(nmout_gray.coeff_stats(3,3)),get_p_html(nmout_gray.coeff_stats(3,4)),...
    (exp(nmout_gray.coeff_stats(4,1))-1)*100,exp(nmout_gray.coeff_stats(4,2)),exp(nmout_gray.coeff_stats(4,3)),get_p_html(nmout_gray.coeff_stats(4,4)));

%% initialize figure
figure
set(gcf,'position',[10 10 1000 1200])
tiledlayout(3,2,'tilespacing','tight','padding','tight')


%% LOO
pt_stats = nmout.pt_stats;
pt_specific = nmout.pt_specific;
aucs = pt_stats(:,8);
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

if 0
    % compare nsoz to n predicted
    table(pred,cellfun(@sum,soz))
end

%% PLot ROC
Y = nmout.Y;
X = nmout.X;
ym = nanmean(Y,1);
std_Y = nanstd(Y,[],1);
ste_Y = std_Y./sqrt(size(Y,1));
ly = ym-std_Y;
uy = ym+std_Y;
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
fprintf(fid,['<p>We next assessed how well this model would '...
    'localize the SOZ in new patients. The mean (standard deviation) AUC of '...
    'an ROC curve representing classifier performance on the leave-one-out test patient was %1.2f (%1.2f) (Fig 5A). '...
    'Results were similar for a model analyzing electrodes in gray matter only '...
    '(mean (SD) AUC =  %1.2f (%1.2f)). '],...
    nanmean(aucs,1),nanstd(aucs,[],1),nanmean(nmout_gray.pt_stats(:,8),1),nanstd(nmout_gray.pt_stats(:,8),[],1));

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
conf_text = {'True negatives','False positives';'False negatives','True positives'};
for ic = 1:2
    for jc = 1:2
        text(ic,jc,sprintf('%s:\n%d',conf_text{jc,ic},all_mat(jc,ic)),'horizontalalignment','center','fontsize',20,'fontweight','bold')
    end
end
title(sprintf(['All patients in aggregate\nAccuracy: %1.1f%%,'...
     ' PPV: %1.1f%%, NPV: %1.1f%%'],...
all_acc*100,...
all_ppv*100,all_npv*100))
set(gca,'fontsize',15)

%% Make sure ppv what i expect
assert(abs(all_ppv-mat(2,2)/(mat(1,2)+mat(2,2)))<0.01)

%% Text
fprintf(fid,['To provide an example of real-world classifier performance, we selected patient-specific thresholds, '...
    'each chosen such that the number of SOZ electrodes predicted by the model equaled the number of '...
    'clinician-defined SOZ electrodes (which varied across patients). Aggregating all patient data, the '...
    'resulting confusion matrix had an accuracy of %1.1f%%, a positive predictive value (PPV) of '...
    '%1.1f%%, and a negative predictive value (NPV) of %1.1f%% (Fig 5B). '],...
    all_acc*100,...
    all_ppv*100,all_npv*100);



%% PPVs and NPVs for each patient
% Sort lowest to highest AUC
nexttile
[~,I] = sort(aucs);
plot(1:length(aucs),sort(aucs),'ko','linewidth',2)
%{
hold on
plot(1:length(ppv),ppv(I),'o','color',[0, 0.4470, 0.7410],'linewidth',2)
hold on
plot(1:length(ppv),npv(I),'o','linewidth',2)
%}
%plot(1:length(aucs),sort(aucs),'ko','linewidth',2)
%legend({'AUC'},'fontsize',15,'location','southeast')
title('Individual patient accuracies')
%title('Individual patient PPV and NPV')
set(gca,'fontsize',15)
xlabel('Patient')
%ylabel('AUC')
ylabel('AUC')

fprintf(fid,'Model performance was highly variable across patients (Fig 5C).');

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
title('AUC by implant strategy')
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
    'classification AUC '...
    '(median %1.2f) to those with grid/strip/depth implantations'...
    ' (median %1.2f) (Mann-Whitney test: <i>U</i>'...
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
title('AUC by seizure localization')
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
    'classification AUC '...
    '(median %1.2f) relative to those with extra-temporal'...
    ' lobe epilepsy (median %1.2f) (Mann-Whitney test: <i>U</i>'...
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
time_aucs = nmout.time_aucs;
sp = nan(2,1);
for iws = 1:2
    data = squeeze(time_aucs(iws,:,:,:));
    
    %% For each duration, do a linear mixed effects model to estimate standard error
    % initialize mean and ste
    mean_ste = nan(size(data,1),2);
    
    % loop over durations
    for id = 1:size(data,1)
        
        curr_aucs = squeeze(data(id,:,:));
        
        [mean_estimate,ste_estimate] = model_clustered_effect(curr_aucs);
        mean_ste(id,:) = [mean_estimate,ste_estimate];
        
    end
    sp(iws) = errorbar(duration_ticks,mean_ste(:,1),mean_ste(:,2),'-o','linewidth',2);
    %flattened_data = reshape(data,size(data,1),size(data,2)*size(data,3));
    %mean_aucs = nanmean(flattened_data,2);
    %std_aucs = nanstd(flattened_data,[],2);
    %ste_aucs = std_aucs/sqrt(size(flattened_data,2));
    %sp(iws) = plot(duration_ticks,curr_aucs,'-o','linewidth',2);
    %sp(iws) = errorbar(duration_ticks,mean_aucs,ste_aucs,'-o','linewidth',2);
    hold on
end
legend({'Wake','Sleep'},'fontsize',15,'location','southeast')
xticks(duration_ticks)
xticklabels(duration_ticklabels)
ylabel('Mean (standard error) AUC')
xlabel('Duration examined (minutes)')
title('Accuracy by duration')
set(gca,'fontsize',15);

% Text
%
sleep_data_10 = squeeze(time_aucs(2,3,:,:));
sleep_data_full = squeeze(time_aucs(2,end,:,:));
[mean_estimate_10,ste_estimate_10] = model_clustered_effect(sleep_data_10);
[mean_estimate_full,ste_estimate_full] = model_clustered_effect(sleep_data_full);
fprintf(fid,['<p>We next studied what duration of interictal data was needed '...
    'to localize the SOZ. The average AUC '...
    'was higher at each duration for models trained on sleep spike data'...
    ' than those trained on wake spike data. Also, near-optimal model performance was '...
    'achieved with only 10 minutes of interictal data (mean (standard error) AUC for model trained on '...
    'full duration of sleep data = %1.2f (%1.3f) versus ten minutes of sleep '...
    'data = %1.2f (%1.3f)) (Fig 5F).</p>'],...
    mean_estimate_full,ste_estimate_full,...
    mean_estimate_10,ste_estimate_10);
%}

%% Annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.9 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.55 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.55 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.23 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.23 0.1 0.1],'String','F','fontsize',25,'linestyle','none')

print(gcf,[out_folder,'Fig5'],'-depsc')
close all

fclose(fid);
fclose(sfid);

end
