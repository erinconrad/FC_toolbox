function out = individual_run_coherence(file_name,times,show_data,show_montage,name)



%% Get file locs
locations = fc_toolbox_locs;

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
fs = data.fs;

%% Cleaned labels
clean_labels = decompose_labels(chLabels,name);

%% Find non-intracranial chs
non_intracranial = find_non_intracranial(clean_labels);
which_chs = find(~non_intracranial); % channels to do analysis on

%% Reject bad channels
[bad,details] = identify_bad_chs(values,which_chs,chLabels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

%% CAR montage
[car_values,car_labels] = car_montage(values,which_chs,clean_labels);

is_run_car = ismember((1:length(clean_labels))',which_chs);

%% Calculate network
% Loop over montages
for im = 2
   
    
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
    
    %% Get coherence
    coh = coherence_calc(run_values,fs);
    
    %% Get bandpower
    bp = bp_calc(run_values,fs,skip);
    
    % save
    out.montage(im).name = montage;
    out.montage(im).labels = curr_labels;
    out.montage(im).is_run = is_run;
    out.montage(im).bp = bp;
    out.montage(im).coh = coh;
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
    simple_plot(tout,out,ex_chs,show_montage,out.montage(show_montage).spikes,...
        only_run,skip)
    %pause
    %close(gcf)
    clear tout
end


end