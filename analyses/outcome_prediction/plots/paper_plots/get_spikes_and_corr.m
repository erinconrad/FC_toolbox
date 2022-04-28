function out = get_spikes_and_corr

tw = 2;
step = 0.1;
file_name = 'HUP212_phaseII';
times = [100003 100009];
window = diff(times);
which_net = 'pc';

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

%% Get elec locs
box_path = locations.box_folder;
elec_path = [box_path,'CNT Implant Reconstructions/'];
temp_out = return_elec_locs(pt_name,elec_path);
loc_labels = temp_out.elec_names;
locs = temp_out.locs;
anatomy = temp_out.anatomy;
clear temp_out
clean_loc_labels = decompose_labels(loc_labels,pt_name);

%% Get ieeg data
data = download_ieeg_data(file_name,login_name,pwfile,times,1); % 1 means get lots of data
chLabels = data.chLabels;
values = data.values;
raw_values = values;
fs = data.fs;

%% Cleaned labels
clean_labels = decompose_labels(chLabels,pt_name);

%% Reconcile locs and anatomy with ieeg labels
[ieeg_locs,ieeg_anatomy] = reconcile_locs_ieeg(clean_labels,clean_loc_labels,locs,anatomy);

% I don't need these anymore and they might confuse me. I want everything
% to be referenced to the ieeg channels
clear anatomy
clear locs
clear loc_labels
clear clean_loc_labels

locs = ieeg_locs;
anatomy = ieeg_anatomy;

%% Find non-intracranial chs
non_intracranial = find_non_intracranial(clean_labels);
which_chs = find(~non_intracranial); % channels to do analysis on

%% Reject bad channels
[bad,details] = identify_bad_chs(values,which_chs,chLabels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

%% CAR montage
[car_values,car_labels] = car_montage(values,which_chs,clean_labels);
is_run_car = ismember((1:length(clean_labels))',which_chs);
values = car_values;
is_run = is_run_car;
curr_locs = locs;
curr_labels = car_labels;
curr_anatomy = anatomy;

 % make non run channels nans
values(:,~is_run) = nan;
skip = find(~is_run);

% filters
values = notch_filter(values,fs);
values = bandpass_filter(values,fs);

gdf = detector_alt(values,fs);
sp_chs = unique(gdf(:,1));

%% Get networks
% Define time windows
iw = round(tw*fs);
istep = round(step*fs);
istart = (1:istep:size(values,1)-iw)';
iend = (iw:istep:size(values,1))';
nsteps = length(istart);
all_corrs = nan(nsteps,1);
mean_indices = (istart+iend)/2;
mean_times = mean_indices/fs;


% Get running time window of corrs
for is = 1:nsteps
    clip = values(istart(is):iend(is),sp_chs);
    pc = corrcoef(clip);
    pc(logical(eye(size(pc)))) = nan;
    all_corrs(is) = nanmean(pc,'all');
end

times = linspace(0,window,size(values,1));
out.times = times;
out.sp_chs = sp_chs;
out.values = values;
out.mean_times = mean_times;
out.all_corrs = all_corrs;
out.window = window;
%{
nexttile
offset = 0;
for isp = 1:length(sp_chs)
    plot(times,values(:,sp_chs(isp))-offset,'k')
    hold on
    if isp < length(sp_chs)
        offset = offset - (min(values(:,sp_chs(isp))) - max(values(:,sp_chs(isp+1))));
    end
end

nexttile
plot(mean_times,all_corrs)
xlim([0 window])
%}
end