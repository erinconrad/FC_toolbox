

function [all_sleep,all_wake] = compare_sw_ad

%{ 
This function calculates an ROC curve representing the ability of the alpha
delta ratio to discriminate sleep vs wake (as defined by Erin's review of
scalp eeg data)
%}


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
%data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Get the alpha delta ratios for the sleep/wake annotations
%summ = match_sw_ad(pt,summ);
summ = sw_ad_erin_designations;


%% Get a point for each patient
all_sleep = []; 
all_wake = [];
for p = 1:length(summ)
    
    if ~isfield(summ(p),'sw')
        continue
    end
    
    sleep_ad = (summ(p).sw.sleep);
    wake_ad = (summ(p).sw.wake);
    
    % get all ads
    ad = summ(p).ad; % mean across electrodes

    sleep_norm = (sleep_ad-nanmedian(ad))./iqr(ad);
    wake_norm = (wake_ad-nanmedian(ad))./iqr(ad);

    
    all_sleep = [all_sleep;sleep_norm];
    all_wake = [all_wake;wake_norm];
end

%% Calculate roc
[roc,auc] = calculate_roc(all_sleep,all_wake,1e3);


%
figure
plot(roc(:,1),roc(:,2),'k','linewidth',2)
hold on
plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','northwest')
print([out_folder,'roc'],'-dpng');
close all

end

%{
%% Plot
figure
for p = 1:npts
    if any(isnan(all_sw(p,:))), continue; end
    plot([1 2],[all_norm(p,1),all_norm(p,2)],'k-');
    hold on
end
xlim([0.5 2.5])

[~,pval] = ttest(all_norm(:,1),all_norm(:,2))
%}

