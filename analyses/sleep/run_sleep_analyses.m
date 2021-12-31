%% Overview
%{
This script takes intermediate datasets (summ.mat, one for each patient)
that contain spike counts and network calculation and does analyses to see
what happens to spikes, etc., with sleep. It output an out.mat into the
data folder to be used for plots and statistical tests.
%}

%% Parameters
rm_cluster = 0; % remove clustered seizures (should be no)
do_avg = 0; % should be no
exc = []; % number of blocks to exclude (should be empty)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
out_folder = [scripts_folder,'analyses/sleep/data/'];

addpath(genpath(scripts_folder))

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end
%}

%% Do circadian analysis
fprintf('\nDoing circadian analysis\n'); 
circ_out = all_pt_psd;

%% Do alpha delta ratio validation
fprintf('\nDoing alpha delta ratio validation\n'); 
[roc,auc,disc,disc_I,swdes] = ad_validation(exc);
roc_out.roc = roc;
roc_out.auc = auc;
roc_out.disc = disc;
roc_out.disc_I = disc_I;
roc_out.swdes = swdes;



%% Do binary ad analyses
fprintf('\nDoing binary AD analyses\n');
bin_out = binary_ad_analyses(disc,exc);

% Stats on amount of wake and sleep
n_sleep_wake = bin_out.n_sleep_wake;
perc_asleep = n_sleep_wake(:,1)./(sum(n_sleep_wake,2))*100;
iqr_sleep = prctile(perc_asleep,[25 75]);
fprintf('\nAcross all patients, a median of %1.1f%% (IQR %1.1f%% - %1.1f%%) of periods were determined to be asleep.\n',...
    nanmedian(perc_asleep),iqr_sleep(1),iqr_sleep(2));


% verified through here
%% Seizure circadian analysis
sz_circ_out = sz_circadian(disc,exc);

%% Do seizure time analyses
fprintf('\nDoing seizure histogram analyses\n');
sz_out = seizure_time_histogram(rm_cluster,do_avg,disc);

%% Put together
out.circ_out = circ_out;
out.roc_out = roc_out;
%out.subnet_out = subnet_out;
%out.sleep_hist_out = sleep_hist_out;
out.sz_circ_out = sz_circ_out;
out.bin_out = bin_out;
out.sz_out = sz_out;
out.out_folder = out_folder;

%% Save out file
save([out_folder,'out.mat'],'out')

%% Do plots
%fprintf('\nDoing other plots\n');
%sleep_plots(out,1)

