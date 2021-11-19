%% To do
% Don't need to skip patients in circadian analysis

%% Parameters
rm_cluster = 1;
do_avg = 0;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
int_folder = [results_folder,'analysis/intermediate/'];
out_folder = [results_folder,'analysis/sleep/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

%% Do circadian analysis
fprintf('\nDoing circadian analysis\n');
circ_out = all_pt_psd;

%% Do alpha delta ratio validation
fprintf('\nDoing alpha delta ratio validation\n');
[roc,auc,disc,disc_I] = ad_validation;
roc_out.roc = roc;
roc_out.auc = auc;
roc_out.disc = disc;
roc_out.disc_I = disc_I;

%% Do histogram analysis
fprintf('\nDoing sleep histogram analysis\n');
sleep_hist_out = sleep_histogram_analysis(rm_cluster,disc);

%% Do binary ad analyses
fprintf('\nDoing binary AD analyses\n');
bin_out = binary_ad_analyses(disc);

% Stats on amount of wake and sleep
n_sleep_wake = bin_out.n_sleep_wake;
perc_asleep = n_sleep_wake(:,1)./(sum(n_sleep_wake,2))*100;
iqr_sleep = prctile(perc_asleep,[25 75]);
fprintf('\nAcross all patients, a median of %1.1f%% (IQR %1.1f%% - %1.1f%%) of periods were determined to be asleep.\n',...
    nanmedian(perc_asleep),iqr_sleep(1),iqr_sleep(2));

%% Do seizure time analyses
fprintf('\nDoing seizure histogram analyses\n');
sz_out = seizure_time_histogram(rm_cluster,do_avg,disc);

%% Put together
out.circ_out = circ_out;
out.roc_out = roc_out;
out.sleep_hist_out = sleep_hist_out;
out.bin_out = bin_out;
out.sz_out = sz_out;
out.out_folder = out_folder;

%% Do plots
fprintf('\nDoing other plots\n');
sleep_plots(out)

%% Save out file
save([out_folder,'out.mat'],'out')