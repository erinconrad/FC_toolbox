function compare_sw_ad(pt,summ)

%% Get the alpha delta ratios for the sleep/wake annotations
summ = match_sw_ad(pt,summ);


%% Get a point for each patient
npts = length(summ);
all_sw = nan(npts,2); %1 = sleep, 2 = wake
all_prctile = nan(npts,2); %1 = sleep, 2 = wake
all_norm = nan(npts,2); 
for p = 1:npts
    sleep_ad = nanmean(summ(p).sw.sleep);
    wake_ad = nanmean(summ(p).sw.wake);
    
    all_sw(p,:) = [sleep_ad wake_ad];

    % get all ads
    ad = nanmean(summ(p).ad,1); % mean across electrodes
    
    % get the percentile
    sleep_prctile = mean(ad <= sleep_ad);
    wake_prctile = mean(ad <= wake_ad);

    sleep_norm = (sleep_ad-nanmean(ad))./nanstd(ad);
    wake_norm = (wake_ad-nanmean(ad))./nanstd(ad);

    all_prctile(p,:) = [sleep_prctile wake_prctile];
    all_norm(p,:) = [sleep_norm,wake_norm];
end



%% Plot
figure
for p = 1:npts
    if any(isnan(all_sw(p,:))), continue; end
    plot([1 2],[all_norm(p,1),all_norm(p,2)],'k-');
    %plot([1 2],[all_sw(p,1),all_sw(p,2)],'k-');
    %plot([1 2],[all_prctile(p,1),all_prctile(p,2)],'k-');
    hold on
end
xlim([0.5 2.5])

%[~,pval] = ttest(all_prctile(:,1),all_prctile(:,2))
%[~,pval] = ttest(all_sw(:,1),all_sw(:,2))
[~,pval] = ttest(all_norm(:,1),all_norm(:,2))

end