%% Plot desired spikes

% which pts
%pts = {'HUP162','HUP157','MP0001','MP0003'};
pts = {'HUP157'};
npts = length(pts);
im = 3; %bipolar

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
edf_path = [results_folder,'edf_out/'];
edf_summ_path = [results_folder,'edf_summ_out_epilepsy_laterality/'];
sleep_stage_path = [results_folder,'edf_out/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

validation_file = [scripts_folder,'spike_detector/Manual validation.xlsx'];
mT = readtable(validation_file,'Sheet','strange_elec_names');


%% Loop over patients
for p = 1:npts

    name = pts{p};

    % Load the edf summary file
    info = load([edf_summ_path,name,'/summ.mat']);
    info = info.out;

    random_spikes_pretty(info.all_spike_times{im},name,info.labels,info.montages{im},edf_path,edf_summ_path,mT)

end
