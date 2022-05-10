function T1 = main_table_fc


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/plots/paper_plots/'];
int_folder = [results_folder,'analysis/intermediate/']; % this includes patients with good and bad spikes
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

names = cell(npts,1);
age_onset = nan(npts,1);
age_implant = nan(npts,1);
sex = cell(npts,1);
nelecs = nan(npts,1);
rate = nan(npts,1);
duration = nan(npts,1);
lat = cell(npts,1);
loc = cell(npts,1);
stereo = nan(npts,1);

%% Loop over files
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Basic demographics
    names{p} = summ.name;
    clinical = summ.clinical;
    age_onset(p) = clinical.age_onset;
    sex{p} = clinical.sex;
    age_implant(p) = clinical.age_implant;
    stereo(p) = clinical.stereo;
    
    %% Electrodes
    labels = summ.labels;
    ekg = find_non_intracranial(labels);
    labels = labels(~ekg);
    
    nelecs(p) = length(labels);

    %% Duration
    duration(p) = summ.times(end)/3600/24;
    
    %% SOZ localization
    %[curr_loc,curr_lat] = seizure_localization_parser(summ.soz.loc,summ.soz.lat);
    [curr_loc,curr_lat] = bin_manual_locs(summ.soz.loc,summ.soz.lat);
    lat{p} = curr_lat;
    loc{p} = curr_loc;
    
    
end

T1 = table(names,sex,age_onset,age_implant,nelecs,...
    duration,rate,loc,lat);

%% Turn into summary stats

% Get data
nfemale = sum(cellfun(@(x) strcmp(x,'Female'),sex));
median_range_age_onset = [nanmedian(age_onset),min(age_onset),max(age_onset)];
median_range_age_implant = [nanmedian(age_implant),min(age_implant),max(age_implant)];
median_range_nelecs = [nanmedian(nelecs),min(nelecs),max(nelecs)];
median_range_duration = [nanmedian(duration),min(duration),max(duration)];
n_left = sum(cellfun(@(x) strcmp(x,'left'),lat));
n_right = sum(cellfun(@(x) strcmp(x,'right'),lat));
n_bilateral = sum(cellfun(@(x) strcmp(x,'bilateral'),lat));
n_mesial_temporal = sum(cellfun(@(x) strcmp(x,'mesial temporal'),loc));
n_temporal_neocortical = sum(cellfun(@(x) strcmp(x,'temporal neocortical'),loc));
n_other_cortex = sum(cellfun(@(x) strcmp(x,'other cortex'),loc));

% Turn to table
total_str = {'Total: N',sprintf('%d',npts)};
female_str = {'Female: N (%)',sprintf('%d (%1.1f%%)',nfemale,nfemale/npts*100)};
age_onset_str = {'Age at onset in years: median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_age_onset(1),median_range_age_onset(2),median_range_age_onset(3))};
age_implant_str = {'Age at implant in years: median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_age_implant(1),median_range_age_implant(2),median_range_age_implant(3))};
nelecs_str = {'Number of electrodes: median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_nelecs(1),median_range_nelecs(2),median_range_nelecs(3))};
implant_str = {'Implant type',''};
n_gs_str = {'Grids/strips/depths: N (%)',sprintf('%d (%1.1f%%)',sum(stereo==0),sum(stereo==0)/length(stereo)*100)};
n_stereo_str = {'Stereo-EEG: N (%)',sprintf('%d (%1.1f%%)',sum(stereo==1),sum(stereo==1)/length(stereo)*100)};
duration_str = {'Intracranial recording duration in days: median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_duration(1),median_range_duration(2),median_range_duration(3))};
lat_str = {'Seizure laterality',''};
left_str = {'Left: N (%)',sprintf('%d (%1.1f%%)',n_left,n_left/npts*100)};
right_str = {'Right: N (%)',sprintf('%d (%1.1f%%)',n_right,n_right/npts*100)};
bilateral_str = {'Bilateral: N (%)',sprintf('%d (%1.1f%%)',n_bilateral,n_bilateral/npts*100)};
loc_str = {'Seizure localization',''};
mesial_temporal_str = {'Mesial temporal: N (%)',sprintf('%d (%1.1f%%)',n_mesial_temporal,n_mesial_temporal/npts*100)};
temporal_neocortical_str = {'Temporal neocortical: N (%)',sprintf('%d (%1.1f%%)',n_temporal_neocortical,n_temporal_neocortical/npts*100)};
other_cortex_str = {'Other cortex: N (%)',sprintf('%d (%1.1f%%)',n_other_cortex,n_other_cortex/npts*100)};

%{
Rows:
Total N
Female/Male
Age at onset
Age at implant
Seizure laterality
Seizure localization - MT, T, Other cortex, multifocal
N g/s/d vs stereo
Median (range) # electrodes
%}


all = [total_str;...
    female_str;...
    age_onset_str;...
    age_implant_str;...
    lat_str;...
    left_str;...
    right_str;...
    bilateral_str;...
    loc_str;...
    mesial_temporal_str;...
    temporal_neocortical_str;
    other_cortex_str;...
    nelecs_str;...
    implant_str;...
    n_gs_str;...
    n_stereo_str;...
    duration_str];

%{
all = [total_str;...
    female_str;...
    age_onset_str;...
    age_implant_str;...
    nelecs_str;...
    implant_str;...
    n_gs_str;...
    n_stereo_str;...
    duration_str;...
    rate_str;...
    lat_str;...
    left_str;...
    right_str;...
    bilateral_str;...
    loc_str;...
    temporal_str;...
    other_str];
%}

T2 = cell2table(all);
writetable(T2,[out_folder,'Table1.csv']);

end



