%% Parameters
file_name = 'HUP78_phaseII-Annotations';%'HUP212_phaseII';
times = [3.0480e+05-0 3.0480e+05+7];

%% get pt name
C = strsplit(file_name,'_');
pt_name = C{1};

%% Get file locs
locations = fc_toolbox_locs;

% output folder
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'tests/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Get ieeg data
data = download_ieeg_data(file_name,login_name,pwfile,times,1); % 1 means get lots of data

chLabels = data.chLabels;
values = data.values;
raw_values = values;
fs = data.fs;

%% Cleaned labels
clean_labels = decompose_labels(chLabels,pt_name);

%% Find non-intracranial chs
non_intracranial = find_non_intracranial(clean_labels);
which_chs = find(~non_intracranial); % channels to do analysis on

%% Reject bad channels
[bad,details] = identify_bad_chs(values,which_chs,chLabels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

%% CAR montage
[car_values,car_labels] = car_montage(values,which_chs,clean_labels);
values = car_values;
labels = car_labels;
is_run = ismember((1:length(clean_labels))',which_chs);

%% filters
values = notch_filter(values,fs);
values = bandpass_filter(values,fs);

%%
run_values = values;
run_values(:,~is_run) = nan;

%% Spike detectors
gdf = detector_new_timing(run_values,fs);

%%
show_chs = ismember(labels,{'LG56-CAR','LG55-CAR','LG62-CAR','LG63-CAR','LG54-CAR','LG46-CAR','LG47-CAR'});

%% Reorder by spike time
new_gdf = gdf(ismember(gdf(:,1),find(show_chs)),:);
[~,I] = sort(new_gdf(:,2),'ascend');

show_eeg_and_spikes(values,labels,gdf,fs,new_gdf(I,1))