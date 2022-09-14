function validate_raw_ad

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% get my manual sleep/wake designations
swdes = sw_ad_erin_designations;

npts_val = length(swdes);
ad = nan(npts_val,2); %1 = sleep, 2 = wake
all_wake = [];
all_sleep = [];

all_vals = [];
all_labels = {};
for j = 1:npts_val
    if isempty(swdes(j).sw), continue; end
    sleep_ad = swdes(j).sw.sleep;
    wake_ad = swdes(j).sw.wake;

    
    ad(j,:) = [nanmean(sleep_ad),nanmean(wake_ad)]; 
    
    % note that for generating the roc curve, I combine ALL values (rather
    % than take an average value per patient). This gets me more datapoints
    % but it would give some preference to the patients for whom I could
    % successfully discern more sleep/wake periods (which isn't
    % unreasonable).
    all_wake = [all_wake;wake_ad];
    all_sleep = [all_sleep;sleep_ad];
    
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