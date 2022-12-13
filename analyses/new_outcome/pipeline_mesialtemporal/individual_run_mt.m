function out = individual_run_mt(file_path)

%% Parameters
show_data = 0;
tw = 2;

%% Load the edf file
C = strsplit(file_path,'/');
name = C{end-1};
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
[bad,details] = identify_bad_chs(curr_values,which_chs,allowed_labels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

%% CAR reference
[car_values,car_labels] = car_montage(curr_values,which_chs,allowed_labels);
is_run_car = ismember((1:length(car_labels))',which_chs);

%% Machine reference
machine_values = curr_values;
machine_labels = cellfun(@(x) sprintf('%s-Ref',x),allowed_labels,'uniformoutput',false);
is_run_machine = ismember((1:length(car_labels))',which_chs);

%% Bipolar reference
[bipolar_values,~,bipolar_labels,chs_in_bipolar] = ...
    bipolar_montage_fc(curr_values,allowed_labels,[],[],name);
bad_bipolar = any(ismember(chs_in_bipolar,bad),2);
empty = cellfun(@(x) strcmp(x,'-'),bipolar_labels);
which_chs_bipolar = 1:size(chs_in_bipolar,1);
which_chs_bipolar(bad_bipolar|empty) = [];
is_run_bipolar = ismember((1:length(allowed_labels))',which_chs_bipolar);

%% Calculate network
% Loop over montages
for im = 1:3 
   
    if im == 1
        montage = 'machine';
        values = machine_values;
        curr_labels = machine_labels;
        is_run = is_run_machine;
    elseif im == 2
        montage = 'car';
        values = car_values;
        is_run = is_run_car;
        curr_labels = car_labels;
    elseif im == 3
        montage = 'bipolar';
        values = bipolar_values;
        is_run = is_run_bipolar;
        curr_labels = bipolar_labels;
    end

    % filters
    values = notch_filter(values,fs);
    values = bandpass_filter(values,fs);
    
    % make non run channels nans
    run_values = values;
    run_values(:,~is_run) = nan;
    skip = find(~is_run);
    
    % PC
    pc =  wrap_or_unwrap_adjacency_fc_toolbox(pc_vector_calc(run_values,fs,tw));
    
    % Spikes
    gdf = detector_new_timing(run_values,fs);
    
    % Get alpha delta ratio
    ad_rat = calc_ad(run_values,fs);
    
    % Get bandpower
    bp = bp_calc(run_values,fs,[]);
    
    % Get coherence    
    coh = faster_coherence_calc(run_values,fs);

    % PLV
    plv = plv_calc(run_values,fs);
    
    out.montage(im).name = montage;
    out.montage(im).bp = bp;
    out.montage(im).pc = pc;
    out.montage(im).plv = plv;
    out.montage(im).coh = coh;
    out.montage(im).gdf = gdf;
    out.montage(im).ad = ad_rat;
    out.montage(im).skip = skip;
    out.montage(im).is_run = is_run;
    out.montage(im).labels = curr_labels;
    out.clean_labels = allowed_labels;
    out.fs = fs;
    out.times = [curr_times(1) curr_times(end)];
    out.idx = [rand_start rand_end];
    out.file_path = file_path;
    out.name = name;

    if show_data
        tout.montage(im).values = values;
        tout.montage(im).name = montage;
    end

end


 if show_data
    ex_chs = [];
    only_run = 0;
    show_montage = 3;
    simple_plot(tout,out,ex_chs,show_montage,out.montage(show_montage).gdf,...
        only_run,skip)
    %pause
    %close(gcf)
    clear tout

    figure; tiledlayout(1,3);
    for i = 1:3
        nexttile
        turn_nans_gray(out.montage(i).pc)
        title(out.montage(i).name)
        xticks(1:length(out.montage(i).labels))
        xticklabels(out.montage(i).labels)
        yticks(1:length(out.montage(i).labels))
        yticklabels(out.montage(i).labels)
    end

    figure
    tiledlayout(1,3)
    for i = 1:3
        nexttile
        turn_nans_gray(out.montage(i).plv(:,:,5))
        title(out.montage(i).name)
        xticks(1:length(out.montage(i).labels))
        xticklabels(out.montage(i).labels)
        yticks(1:length(out.montage(i).labels))
        yticklabels(out.montage(i).labels)
    end
end


end