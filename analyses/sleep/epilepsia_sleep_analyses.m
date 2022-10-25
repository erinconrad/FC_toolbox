%% Overview
%{
This script takes intermediate datasets (summ.mat, one for each patient)
that contain spike counts and network calculation and does analyses to see
what happens to spikes, etc., with sleep. It output an out.mat into the
data folder to be used for plots and statistical tests.
%}

%% Parameters to change
doing_from_github = 1; % change to 1 if doing from github

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
out_folder = [scripts_folder,'analyses/sleep/data/'];
addpath(genpath(scripts_folder))

%% Do circadian analysis
fprintf('\nDoing circadian analysis\n'); 
circ_out = all_pt_psd;

%% Do alpha delta ratio validation and get ADR cutoff for wake/sleep
fprintf('\nDoing alpha delta ratio validation\n'); 
roc_out = ad_validation;
disc = roc_out.disc;

%% Do SleepSEEG analyses
fprintf('\nDoing SleepSEEG analyses\n')
seeg_out = basic_seeg_analyses;
seeg_ad_out = sleep_seeg_ad(0);

%% Time-varying analysis
fprintf('\nDoing time-varying analysis\n');
time_out = time_varying_spikes(disc);

%% Do binary ad analyses
fprintf('\nDoing binary AD analyses\n');
bin_out = binary_ad_analyses(disc);

% Stats on amount of wake and sleep
n_sleep_wake = bin_out.n_sleep_wake;
perc_asleep = n_sleep_wake(:,1)./(sum(n_sleep_wake,2))*100;
iqr_sleep = prctile(perc_asleep,[25 75]);
fprintf('\nAcross all patients, a median of %1.1f%% (IQR %1.1f%% - %1.1f%%) of periods were determined to be asleep.\n',...
    nanmedian(perc_asleep),iqr_sleep(1),iqr_sleep(2));


%% Seizure circadian analysis
sz_circ_out = sz_circadian(disc);

%% Do seizure time analyses
fprintf('\nDoing seizure histogram analyses\n');
sz_out = seizure_time_histogram(disc);

%% Put together
out.circ_out = circ_out;
out.roc_out = roc_out;
out.sz_circ_out = sz_circ_out;
out.bin_out = bin_out;
out.sz_out = sz_out;
out.time_out = time_out;
out.out_folder = out_folder;
out.seeg_out = seeg_out;
out.seeg_ad_out = seeg_ad_out;

save([out_folder,'out.mat'],'out')

%% Do sleep models
% Skip this if I am only doing the periictal analysis
if exist('just_for_periictal','var') ~= 0 && just_for_periictal == 1
    return
end
fprintf('\nDoing sleep models\n');
fprintf('\nDoing all electrodes model\n');
model_out = sleep_models(0);

fprintf('\nDoing gray matter electrodes model\n');
model_out_gray = sleep_models(1);



out.model_out = model_out;
out.model_out_gray = model_out_gray;


%% Save out file
save([out_folder,'out.mat'],'out')

%% Do plots
if 0
fprintf('\nDoing plots\n');
if ~doing_from_github
    epilepsia_figure1 % this needs a raw data file to run (which is large)
end

epilepsia_figure2

epilepsia_figure3

epilepsia_figure4

epilepsia_figure5

epilepsia_supplemental_fig1

epilepsia_supplemental_fig2

if ~doing_from_github
    epilepsia_table1 % this needs a raw data file to run (which is large)
end
end

