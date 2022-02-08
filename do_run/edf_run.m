function edf_run

%% Parameters
data_folder = '/Users/erinconrad/Desktop/research/other people work/jeremy spikes/';
data_file = 'P004_Seizure_Natus.edf';
name = 'P004';
chan_names = ["Fp1" "Fp2" "F3" "F4" "F7" "F8" "C3" "C4" "Fz" "Cz" "T7" ...
    "T8" "P7" "P8" "O1" "O2"];
which_net = 'pc';
show_data = 1;
tw = 2;

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% ieeg stuff
%{
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
%}

%% Get ieeg data
%{
data = download_ieeg_data(file_name,login_name,pwfile,times,1); % 1 means get lots of data
chLabels = data.chLabels;
values = data.values;
fs = data.fs;
%}

%% Get data
edf_data = edfread([data_folder,data_file],'SelectedSignals',chan_names);
%edf_data = edfread(filename);
info = edfinfo([data_folder,data_file]);
fs = info.NumSamples/seconds(info.DataRecordDuration);
fs = fs(1);
chLabels = cellstr(chan_names);
which_chs = 1:length(chLabels);

num_channels = 16;

data = [];

for j = 1:num_channels
    channels_data = [];
    channel_data = edf_data.(j);
    [len ~] = size(channel_data);
    for i = 1:len
        epoch_samples = channel_data{i};
        channels_data = cat(1, channels_data,epoch_samples);
    end
    data = cat(2,data,channels_data);
end
values = data;

%% Cleaned labels
clean_labels = decompose_labels(chLabels,name);


%% Reject bad channels
[bad,details] = identify_bad_chs(values,which_chs,chLabels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

%% CAR montage
[car_values,car_labels] = car_montage(values,which_chs,clean_labels);

%% Bipolar montage
[bipolar_values,~,bipolar_labels,chs_in_bipolar,which_chs_bipolar] = ...
    bipolar_montage(values,chLabels,which_chs,[],[],name);

is_run_car = ismember((1:length(clean_labels))',which_chs);
is_run_bipolar = ismember((1:length(clean_labels))',which_chs_bipolar);

%% Calculate network
% Loop over montages
for im = 1:2
   
    
    if im == 1
        montage = 'bipolar';
    elseif im == 2
        montage = 'car';
    end


    switch montage
        case 'bipolar'
            values = bipolar_values;
            is_run = is_run_bipolar;
            curr_labels = bipolar_labels;
        case 'car'
            values = car_values;
            is_run = is_run_car;
            curr_labels = car_labels;
    end
    
    % filters
    values = notch_filter(values,fs);
    values = bandpass_filter(values,fs);
    
    % make non run channels nans
    run_values = values;
    run_values(:,~is_run) = nan;
    skip = find(~is_run);
    
    % Choose network
    switch which_net
        case 'pc'
            curr_net = pc_vector_calc(run_values,fs,tw);
        
    end
    
    % Get spikes
    gdf = detector_for_jeremy(run_values,fs);
    fprintf('\nDetected %d spikes\n',size(gdf,1));
    
    if ~isempty(gdf)
        % Detect bad spikes
       % bad_gdf = detector_alt(run_values,fs);
        
        
    end
    
    % Get alpha delta ratio
    ad_rat = calc_ad(run_values,fs);
    
    % save
    out.montage(im).name = montage;
    out.montage(im).net.name = which_net;
    out.montage(im).net.data = curr_net;
    out.montage(im).labels = curr_labels;
    out.montage(im).is_run = is_run;
    out.montage(im).spikes = gdf;
    out.montage(im).ad = ad_rat;
    out.fs = fs;
    out.montage(im).skip = skip;
    out.clean_labels = clean_labels;
    
    if show_data
        tout.montage(im).values = values;
        tout.montage(im).name = montage;
    end
    
    
        
        

end

if show_data
    ex_chs = [];
    only_run = 0;
    show_montage = 1;
    simple_plot(tout,out,ex_chs,show_montage,out.montage(show_montage).spikes,...
        only_run,skip)
    %pause
    %close(gcf)
    clear tout
end


end