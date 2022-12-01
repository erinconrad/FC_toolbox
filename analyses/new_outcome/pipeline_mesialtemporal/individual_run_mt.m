function out = individual_run_mt(file_path)

%% Parameters
tw = 2;

%% Load the edf file
info = edfinfo(file_path);

%% Get basic info
samples_per_time = info.NumSamples(1,1);
num_samples = samples_per_time*info.NumDataRecords;
fs = round(info.NumSamples(1,1)/seconds(info.DataRecordDuration));
labels = cellstr(info.SignalLabels);

%% get allowable electrodes
allowable_labels = get_allowable_elecs;

%% Find labels that match allowable electrodes
allowed = ismember(labels,allowable_labels);
allowed_labels = labels(allowed);
nallowed = sum(allowed);

%% Skip if not left and right
if sum(contains(allowed_labels,'L'))==0 || sum(contains(allowed_labels,'R'))==0
    out = [];
    %return
end

%% Initialize values
values = nan(num_samples,nallowed);

% Separately call edfread for each signal
for is = 1:nallowed
    curr_signal = allowed_labels{is};
    
    % Get timetable for that signal
    T = edfread(file_path,'SelectedSignals',curr_signal);
      
    % Get variable name (may be changed once put into table)
    Var1 = T.Properties.VariableNames{1};
    
    %% Convert the time table to an array
    temp_values = nan(num_samples,1);

    % Loop over segments
    for s = 1:size(T,1)
        seg = T.(Var1){s};

        % Where are we in the temp values
        start_idx = (s-1)*samples_per_time+1;
        end_idx = s*samples_per_time;

        % Fill up values
        temp_values(start_idx:end_idx) = seg;
    end
    
    %% Fill up values
    values(:,is) = temp_values;
    
end

%% Get times
nsamples = size(T,1)*samples_per_time;
times = linspace(0,num_samples/fs,nsamples);

%% Take a random one minute segment
max_start = length(times) - fs*60-1; % must start 60 seconds before the end
rand_start = randi(round(max_start));
rand_end = rand_start + round(fs*60);

%% Narrow values down to this
curr_times = times(rand_start:rand_end);
curr_values = values(rand_start:rand_end,:);

%% Reject bad channels
which_chs = 1:nallowed;
bad = identify_bad_chs(curr_values,which_chs,allowed_labels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

%% CAR montage
[car_values,car_labels] = car_montage(curr_values,which_chs,allowed_labels);
values = car_values;
curr_labels = car_labels;
is_run = ismember((1:length(curr_labels))',which_chs);

% filters
values = notch_filter(values,fs);
values = bandpass_filter(values,fs);

% make non run channels nans
run_values = values;
run_values(:,~is_run) = nan;
skip = find(~is_run);

%% PC
pc =  wrap_or_unwrap_adjacency_fc_toolbox(pc_vector_calc(run_values,fs,tw));

%% Spikes
gdf = detector_new_timing(run_values,fs);

%% Get alpha delta ratio
ad_rat = calc_ad(run_values,fs);

%% Get bandpower
bp = bp_calc(run_values,fs,[]);

%% Get coherence    
coh = faster_coherence_calc(run_values,fs);

out.bp = bp;
out.pc = pc;
out.coh = coh;
out.gdf = gdf;
out.pc = pc;
out.ad = ad_rat;
out.fs = fs;
out.skip = skip;
out.clean_labels = allowed_labels;
out.montage_labels = curr_labels;
out.times = [curr_times(1) curr_times(end)];
out.idx = [rand_start rand_end];
out.file_path = file_path;

end