function roc_out = ad_validation

%{
This is the main validation function to test how well the alpha delta ratio
distingiushes wake from sleep. I specifically use a normalized alpha delta
ratio, where I take the ad and subtract the median and divide by the iqr.
%}

%% Parameter not to change
exc = [];

% get my manual sleep/wake designations
swdes = sw_ad_erin_designations;

npts_val = length(swdes);
ad_norm = nan(npts_val,2); %1 = sleep, 2 = wake
all_wake = [];
all_sleep = [];

all_vals = [];
all_labels = {};
for j = 1:npts_val
    if isempty(swdes(j).sw), continue; end
    sleep_ad = swdes(j).sw.sleep;
    wake_ad = swdes(j).sw.wake;
    ad_val = swdes(j).ad;

    %{
    sleep_norm = (sleep_ad-nanmedian(ad_val))./iqr(ad_val);
    wake_norm = (wake_ad-nanmedian(ad_val))./iqr(ad_val);
    %}
    sleep_norm = norm_exc(sleep_ad,ad_val,exc); % If exc is [], which it should be, this is the same as the above commented out thing
    wake_norm = norm_exc(wake_ad,ad_val,exc);
    
    
    ad_norm(j,:) = [nanmean(sleep_norm),nanmean(wake_norm)]; 
    
    % note that for generating the roc curve, I combine ALL values (rather
    % than take an average value per patient). This gets me more datapoints
    % but it would give some preference to the patients for whom I could
    % successfully discern more sleep/wake periods (which isn't
    % unreasonable).
    all_wake = [all_wake;wake_norm];
    all_sleep = [all_sleep;sleep_norm];
    
    all_vals = [all_vals;wake_norm;sleep_norm];
    wake_label = cell(length(wake_norm),1);
    wake_label(:) = {'Wake'};
    sleep_label = cell(length(sleep_norm),1);
    sleep_label(:) = {'Sleep'};
    all_labels = [all_labels;wake_label;sleep_label];
end

% Calculate roc
[roc,auc,disc,disc_I] = calculate_roc(all_sleep,all_wake,1e3);

% alternate approach to roc (double checking AUC)
labels = all_labels;
scores = all_vals;
[X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'Wake');


roc_out.roc = roc;
roc_out.auc = auc;
roc_out.disc = disc;
roc_out.disc_I = disc_I;
roc_out.swdes = swdes;
roc_out.alt_auc = AUC;
roc_out.X = X;
roc_out.Y = Y;
roc_out.T = T;
roc_out.OPTROCPT = OPTROCPT;


end