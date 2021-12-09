function [roc,auc,disc,disc_I,swdes] = ad_validation

%{
This is the main validation function to test how well the alpha delta ratio
distingiushes wake from sleep. I specifically use a normalized alpha delta
ratio, where I take the ad and subtract the median and divide by the iqr.
%}

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
    
    % note that for generating the roc curve, I combine ALL values (rather
    % than take an average value per patient). This gets me more datapoints
    % but it would give some preference to the patients for whom I could
    % successfully discern more sleep/wake periods (which isn't
    % unreasonable).
    all_wake = [all_wake;wake_norm];
    all_sleep = [all_sleep;sleep_norm];
end

% Calculate roc
[roc,auc,disc,disc_I] = calculate_roc(all_sleep,all_wake,1e3);

% Get number of sleep and wake

end