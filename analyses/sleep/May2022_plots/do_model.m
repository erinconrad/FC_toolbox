function do_model(from_scratch)
%{
1) It looks like sleep and post-ictal spike rates independently help localize,
but awake and pre-ictal do not. (CONFIRM)
2) You only need about 30 minutes of sleep to get the most bang for your buck (the other analysis)
3) Patients are heterogeneous in how well spikes localize, with TLE
patients doing much better
4) I need to think of some way to conceptualize actual thresholds. Perhaps
I should set thresholds such that same number of SOZ as median number. How
many false positives and false negatives?? If x is the proportion of
electrodes that I want to be designated SOZ electrodes, that means I want
(TP + FP)/(TP+TP+FN+FP) =  X
I also have the constraint of the ROC curve
%}

rng(0) % seed rng (for bootstrap portion)

%% Parameters

durations = {1, 5, 10, 30, 60, 60*2,[]};
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

%% Get some clinical stuff
stereo = logical(out.circ_out.stereo);
loc = out.circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');

%% Do the stats
if from_scratch
    % Bootstrap stats
    fprintf('\nDoing bootstrap to get model stats\n');
    coeff_stats = sleep_model_bootstrap_stats;

    % LOO
    fprintf('\nDoing LOO analysis to get individual patient performance\n');
    [pt_stats,X,Y,pt_specific] = sleep_loo;

    % Get AUC for sleep and wake as a function of duration
    fprintf('\nDoing duration analysis\n');
    time_aucs = sleep_duration(durations);

    nmout.coeff_stats = coeff_stats;
    nmout.pt_stats = pt_stats;
    nmout.X = X;
    nmout.Y = Y;
    nmout.pt_specific = pt_specific;
    nmout.time_aucs = time_aucs;

    save([time_roc_folder,'new_model_out.mat'],'nmout');
    
else
    nmout = load([time_roc_folder,'new_model_out.mat']);
    nmout = nmout.nmout;
    unpack_any_struct(nmout);
end

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
[npv,ppv,pred,totaln] = cellfun(@(x,y,z) individual_threshold_stats(x,y,z,desired_threshold),...
    scores,soz,threshold,'uniformoutput',false);
npv = cell2mat(npv);
ppv = cell2mat(ppv);
pred = cell2mat(pred);
totaln = cell2mat(totaln);
nanmean(ppv)
nanmean(pred./totaln)

ym = nanmean(Y,1);
ly = prctile(Y,2.5,1);
uy = prctile(Y,97.5,1);

% Plot
nexttile
plot(X,ym,'linewidth',2)
hold on
plot([0 1],[0 1],'k--','linewidth',2)
legend(sprintf('Median AUC %1.2f',nanmedian(aucs)),'location','southeast','fontsize',15)
xlabel('False positive rate')
ylabel('True positive rate')
title('SOZ classification accuracy')
set(gca,'fontsize',15)


%% PPVs and NPVs for each patient
% Sort lowest to highest AUC
nexttile
[~,I] = sort(ppv);
plot(1:length(ppv),ppv(I),'o','linewidth',2)
hold on
plot(1:length(ppv),npv(I),'o','linewidth',2)
legend({'PPV','NPV'},'fontsize',15,'location','southeast')
title('Individual patient PPVs and NPVs for example threshold')
set(gca,'fontsize',15)
xlabel('Patient')
ylabel('PPV and NPV')

%% Do between-patient analyses

% Implant type
nexttile
plot(1+randn(sum(stereo),1)*0.05,aucs(stereo),'o','linewidth',2)
hold on
plot(2+randn(sum(~stereo),1)*0.05,aucs(~stereo),'o','linewidth',2)
xticks([1 2])
xticklabels({'Stereo-EEG','Grid/strip/depths'});
ylabel('Classification AUC')
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
nylim = [yl(1) yl(1) + 1.15*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k','linewidth',2)
text(1.5,ytext,get_p_text(ranksum(aucs(stereo),aucs(~stereo))),...
    'fontsize',15,'horizontalalignment','center')
ylim(nylim)
xlim([0.5 2.5])
set(gca,'fontsize',15)

% Localization
nexttile
plot(1+randn(sum(temporal),1)*0.05,aucs(temporal),'o','linewidth',2)
hold on
plot(2+randn(sum(extra),1)*0.05,aucs(extra),'o','linewidth',2)
xticks([1 2])
xticklabels({'Temporal lobe epilepsy','Extra-temporal lobe epilepsy'});
ylabel('Classification AUC')
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
nylim = [yl(1) yl(1) + 1.15*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k','linewidth',2)
text(1.5,ytext,get_p_text(ranksum(aucs(temporal),aucs(extra))),...
    'fontsize',15,'horizontalalignment','center')
ylim(nylim)
xlim([0.5 2.5])
set(gca,'fontsize',15)



% Prep duration ticks
duration_ticks = cell2mat(durations(1:ndurs-1));
new_lim = 1.2*(max(duration_ticks)-0)+min(duration_ticks);
duration_ticks = [duration_ticks,new_lim];
duration_ticklabels = arrayfun(@(x) sprintf('%d',x),duration_ticks,'uniformoutput',false);
duration_ticklabels(end) = {'Full'};

nexttile([1 2])
sp = nan(2,1);
for iws = 1:2
    curr_aucs = squeeze(nanmedian(time_aucs(iws,:,:),3));
    sp(iws) = plot(duration_ticks,curr_aucs,'-o','linewidth',2);
    hold on
end
legend({'Wake','Sleep'},'fontsize',15,'location','southeast')
xticks(duration_ticks)
xticklabels(duration_ticklabels)
ylabel('Median AUC')
xlabel('Duration examined (minutes)')
title('Accuracy by duration')
set(gca,'fontsize',15);

%% Print figure
print(gcf,[time_roc_folder,'new_fig'],'-dpng')

%{
%% Get proportion of electrodes that are soz
nsoz = cellfun(@(x) sum(x==1),all_soz);
nelecs = cellfun(@length,all_soz);
median_proportion = median(nsoz./nelecs);

%% Get some clinical stuff
stereo = logical(out.circ_out.stereo);
loc = out.circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');

%% Do main model both to get estimates of coefficients and AUCs
all_auc = nan(nb,1);
all_coeffs = nan(ncoeffs,nb);
all_roc = cell(nb,1);
for ib = 1:nb
     fprintf('\nib = %d of %d\n',ib,nb)
     mout = updated_classifier_may2022([],1);
     all_auc(ib) = mout.AUC;
     for ic = 1:ncoeffs
            % get the estimate of the model coefficient
            all_coeffs(ic,ib) = mout.model.Coefficients{ic+1,2}; 
         
     end
     all_roc{ib} = [mout.X,mout.Y];
    
end
[X,Y] = unify_roc(all_roc);

%% Calculate thresholds
% Find the point along the ROC curve that would return about 5% positivity
% rate
closest_point = find_point_roc_proportion(X,mean(Y,1),median_proportion);


% Get confusion matrices
cout = confusion_matrix(predicted,actual,do_plot);

%% Determine bootstrap stats for each coefficient
coeff_stats = nan(ncoeffs,4); % mean, lower CI, higher CI, p
for ic = 1:ncoeffs
    tout = bootstrap_ci_and_p(squeeze(all_coeffs(ic,:)));
    coeff_stats(ic,:) = [tout.mean,tout.CI_95,tout.p];
    
end


%% LOO analysis - get distribution of individual patient AUCs
pt_auc = nan(npts,1);
for ip = 1:npts
    fprintf('\npatient = %d of %d\n',ip,npts);
    curr_soz = all_soz{ip};
    if sum(curr_soz) == 0
        pt_auc(ip) = nan;
    else
        mout = updated_classifier_may2022(ip,1);
        pt_auc(ip) = mout.AUC;
    end
end

% Stereo doesn't matter
if 0
    figure
    plot(1+randn(sum(stereo),1)*0.05,pt_auc(stereo),'o','linewidth',2)
    hold on
    plot(2+randn(sum(~stereo),1)*0.05,pt_auc(~stereo),'o','linewidth',2)

end

% TLE vs eTLE DOES matter
if 0
    figure
    plot(1+randn(sum(temporal),1)*0.05,pt_auc(temporal),'o','linewidth',2)
    hold on
    plot(2+randn(sum(extra),1)*0.05,pt_auc(extra),'o','linewidth',2)

end

%}



end