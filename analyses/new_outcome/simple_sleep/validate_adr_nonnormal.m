function validate_adr_nonnormal

% get my manual sleep/wake designations
swdes = sw_ad_erin_designations;

npts_val = length(swdes);

all_vals = [];
all_labels = {};

for j = 1:npts_val
    if isempty(swdes(j).sw), continue; end
    sleep_ad = swdes(j).sw.sleep;
    wake_ad = swdes(j).sw.wake;
     
    all_vals = [all_vals;wake_ad;sleep_ad];
    wake_label = cell(length(wake_ad),1);
    wake_label(:) = {'Wake'};
    sleep_label = cell(length(sleep_ad),1);
    sleep_label(:) = {'Sleep'};
    all_labels = [all_labels;wake_label;sleep_label];
end

% alternate approach to roc (double checking AUC)
labels = all_labels;
scores = all_vals;
[X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'Wake');

roc_out.alt_auc = AUC;
roc_out.X = X;
roc_out.Y = Y;
roc_out.T = T;
roc_out.OPTROCPT = OPTROCPT;
disc = T((X==OPTROCPT(1))&(Y==OPTROCPT(2)));
roc_out.disc = disc;

if 0
    figure
    unpaired_plot(scores(strcmp(labels,'Wake')),scores(strcmp(labels,'Sleep')),{'Wake','Sleep'},'ADR')
    hold on
    plot([1 2],[disc disc],'k','linewidth',2)
end

end