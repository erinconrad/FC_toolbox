function quick_run

%% Parameters
file_name = 'HUP212_phaseII';
times = [583265 583280];
which_net = 'pc';
tw = 2;
out_name = 'test1';

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
clean_loc_labels = decompose_labels(loc_labels);

%% Get ieeg data
data = download_ieeg_data(file_name,login_name,pwfile,times,1); % 1 means get lots of data
chLabels = data.chLabels;
values = data.values;
fs = data.fs;

%% Cleaned labels
clean_labels = decompose_labels(chLabels);

%% Reconcile locs and anatomy with ieeg labels
[ieeg_locs,ieeg_anatomy] = reconcile_locs_ieeg(clean_labels,clean_loc_labels,locs,anatomy);

% I don't need these anymore and they might confuse me. I want everything
% to be referenced to the ieeg channels
clear anatomy
clear locs
clear loc_labels
clear clean_loc_labels


%% Find non-intracranial chs
non_intracranial = find_non_intracranial(clean_labels);
which_chs = find(~non_intracranial); % channels to do analysis on

%% Reject bad channels
[bad,details] = identify_bad_chs(values,which_chs,chLabels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

%% CAR montage
car_values = car_montage(values,which_chs);

%% Bipolar montage
[bipolar_values,clean_labels,bipolar_labels,chs_in_bipolar,which_chs_bipolar,mid_locs,mid_anatomy] = ...
    bipolar_montage(values,chLabels,which_chs,ieeg_locs,ieeg_anatomy);

%% Table of channels
is_run_car = ismember((1:length(clean_labels))',which_chs);
is_run_bipolar = ismember((1:length(clean_labels))',which_chs_bipolar);
T = table(clean_labels,is_run_car,bipolar_labels,is_run_bipolar,mid_anatomy);

%% Calculate network
% Loop over montages
for im = 1:2
    
    if im == 1
        values = bipolar_values;
    elseif im == 2
        values = car_values;
    end
    
    % Choose network
    switch which_net
        case 'pc'
            net = pc_calc(values,fs,tw);
    end
    
    % save
    if im == 1
        bipolar_net = net;
    elseif im == 2
        car_net = net;
    end
    
end

%% Save the networks
out.net.bp = bipolar_net;
out.net.car = car_net;
out.file_name = file_name;
out.fs = fs;
out.chLabels = chLabels;
out.bipolar_labels = bipolar_labels;
out.clean_labels = clean_labels;
out.is_run_bipolar = is_run_bipolar;
out.is_run_car = is_run_car;
out.chs_in_bipolar = chs_in_bipolar;
out.bad = bad;
out.bad_details = details;
out.non_intracranial = non_intracranial;
save([out_folder,out_name],'out');

%% Show the networks
if 1
    plot_nets_and_corr(bipolar_net,car_net,bipolar_labels,clean_labels)
    
end


end