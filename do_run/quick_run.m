function quick_run(file_name,times,which_net,do_save,plot_montage)

%% Parameters
tw = 2;
out_name = 'test1';

if isempty(file_name)
    file_name = 'HUP212_phaseII';
end

if isempty(times)
    times = [100000 100015];
end

if isempty(which_net)
    which_net = 'pc';
end

if isempty(do_save)
    do_save = 0;
end

if isempty(plot_montage)
    plot_montage = 1;
end


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

%% Bipolar montage
[bipolar_values,~,bipolar_labels,chs_in_bipolar,which_chs_bipolar,mid_locs,mid_anatomy] = ...
    bipolar_montage(values,chLabels,which_chs,locs,anatomy,pt_name);

%% Table of channels
is_run_car = ismember((1:length(clean_labels))',which_chs);
is_run_bipolar = ismember((1:length(clean_labels))',which_chs_bipolar);
if 0
    
    T = table(clean_labels,is_run_car,bipolar_labels,is_run_bipolar,mid_anatomy);
end

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
            curr_locs = mid_locs;
            curr_labels = bipolar_labels;
            curr_anatomy = mid_anatomy;
        case 'car'
            values = car_values;
            is_run = is_run_car;
            curr_locs = locs;
            curr_labels = car_labels;
            curr_anatomy = anatomy;
    end
    
    % make non run channels nans
    values(:,~is_run) = nan;
    skip = find(~is_run);
    
    % filters
    values = notch_filter(values,fs);
    values = bandpass_filter(values,fs);
    
    tout.montage(im).values = values;
    tout.montage(im).name = montage;
    
    % Loop over networks to run
    for in = 1:2

        if in == 1
            which_net = 'pc';
        elseif in == 2
            which_net = 'inv_dist';
        end

        % Choose network
        switch which_net
            case 'pc'

                curr_net = pc_vector_calc(values,fs,tw);
                
            case 'inv_dist'
                curr_net = inverse_dist(curr_locs);
        end
        
        % Get spikes
        gdf = detector_alt(values,fs);
        fprintf('\nDetected %d spikes\n',size(gdf,1));

        % save
        out.montage(im).name = montage;
        out.montage(im).net(in).name = which_net;
        out.montage(im).net(in).data = curr_net;
        out.montage(im).labels = curr_labels;
        out.montage(im).locs = curr_locs;
        out.montage(im).anatomy = curr_anatomy;
        out.montage(im).is_run = is_run;
        out.montage(im).spikes = gdf;
        out.montage(im).values = values;
        out.fs = fs;

    end
    
    
end
tout.montage(3).name = 'raw';
tout.montage(3).values = raw_values;
out.montage(3).name = 'raw';
out.montage(3).labels = clean_labels;
out.montage(3).is_run = ones(length(clean_labels),1);

%% Save the networks

out.file_name = file_name;
out.fs = fs;
out.chLabels = chLabels;
out.clean_labels = clean_labels;
out.chs_in_bipolar = chs_in_bipolar;
out.bad = bad;
out.bad_details = details;
out.non_intracranial = non_intracranial;
out.tout = tout;
if do_save == 1
    save([out_folder,out_name],'out');
end

%% Show data
if 1
    %ex_chs = {'LA1','LA2','LA3','LA4'};
    
    ex_chs = [];
    only_run = 0;
    simple_plot(tout,out,ex_chs,plot_montage,out.montage(plot_montage).spikes,...
        only_run,skip)
    
end

%% Show the networks
if 0
    plot_nets_and_corr(out,im)
    
end

end


