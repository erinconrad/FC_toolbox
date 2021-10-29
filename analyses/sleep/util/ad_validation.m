function [roc,auc,disc] = ad_validation

swdes = sw_ad_erin_designations;
npts_val = length(swdes);
ad_norm = nan(npts_val,2); %1 = sleep, 2 = wake
all_wake = [];
all_sleep = [];
for j = 1:npts_val
    if isempty(swdes(j).sw), continue; end
    sleep_ad = swdes(j).sw.sleep;
    wake_ad = swdes(j).sw.wake;
    ad_val = swdes(j).ad;

    sleep_norm = (sleep_ad-nanmedian(ad_val))./iqr(ad_val);
    wake_norm = (wake_ad-nanmedian(ad_val))./iqr(ad_val);
    ad_norm(j,:) = [nanmean(sleep_norm),nanmean(wake_norm)];
    all_wake = [all_wake;wake_norm];
    all_sleep = [all_sleep;sleep_norm];
end

% Calculate roc
[roc,auc,disc] = calculate_roc(all_sleep,all_wake,1e3);



end