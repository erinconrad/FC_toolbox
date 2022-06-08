function epilepsia_table1

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/epilepsia/'];
int_folder = [results_folder,'analysis/backup_intermediate_Feb26_good_spikes/'];
%int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
out_folder1 = [scripts_folder,'analyses/sleep/data/'];

%% Load out file
out = load([out_folder1,'out.mat']);
out = out.out;
n_sleep_wake = out.bin_out.n_sleep_wake;

%% Load validation file
val_T = readtable(['Manual validation.xlsx']);
sw_val_T = readtable(['Manual validation.xlsx'],'Sheet','Validation Sw');

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

names = cell(npts,1);
age_onset = nan(npts,1);
age_implant = nan(npts,1);
sex = cell(npts,1);
nelecs = nan(npts,1);
%any_grids = nan(npts,1);
rate = nan(npts,1);
duration = nan(npts,1);
lat = cell(npts,1);
loc = cell(npts,1);
stereo = nan(npts,1);
ppv = nan(npts,1);
perc_asleep = n_sleep_wake(:,1)./sum(n_sleep_wake,2)*100;
ws_ppv = nan(npts,2);

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
    %any_grids(p) = decide_if_any_grids_or_strips(labels);
    
    %% Spike rate
    spikes = summ.spikes;
    rate(p) = nanmean(spikes,'all');
    
    
    %% Duration
    duration(p) = summ.times(end)/3600/24;
    
    %% SOZ localization
    [curr_loc,curr_lat] = seizure_localization_parser(summ.soz.loc,summ.soz.lat);
    lat{p} = curr_lat;
    loc{p} = curr_loc;
    
     %% Spike PPV
    row = strcmp(val_T.name,names);
    assert(sum(row) == 1)
    ppv(p) = val_T.PPV_car_(row);
    
    %% Sleep wake PPV
    row = strcmp(sw_val_T.name,names);
    assert(sum(row) == 1)
    ws_ppv(p,:) = [sw_val_T.x_Correct_outOf50_Wake(row)/50, sw_val_T.x_Correct_outOf50_Sleep(row)/50];
    
    
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
median_range_rate = [nanmedian(rate),min(rate),max(rate)];
median_range_ppv = [nanmedian(ppv),min(ppv),max(ppv)];
median_range_sleep = [nanmedian(perc_asleep),min(perc_asleep),max(perc_asleep)];
n_left = sum(cellfun(@(x) strcmp(x,'left'),lat));
n_right = sum(cellfun(@(x) strcmp(x,'right'),lat));
n_bilateral = sum(cellfun(@(x) strcmp(x,'bilateral'),lat));
n_temporal = sum(cellfun(@(x) strcmp(x,'temporal'),loc));
n_other = sum(cellfun(@(x) strcmp(x,'other'),loc));

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
rate_str = {'Spike rate (spikes/elecs/min): median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_rate(1),median_range_rate(2),median_range_rate(3))};
lat_str = {'Seizure laterality',''};
left_str = {'Left: N (%)',sprintf('%d (%1.1f%%)',n_left,n_left/npts*100)};
right_str = {'Right: N (%)',sprintf('%d (%1.1f%%)',n_right,n_right/npts*100)};
bilateral_str = {'Bilateral: N (%)',sprintf('%d (%1.1f%%)',n_bilateral,n_bilateral/npts*100)};
loc_str = {'Seizure localization',''};
temporal_str = {'Temporal: N (%)',sprintf('%d (%1.1f%%)',n_temporal,n_temporal/npts*100)};
other_str = {'Extra-temporal: N (%)',sprintf('%d (%1.1f%%)',n_other,n_other/npts*100)};
ppv_str = {'Spike detector positive predictive value: median (range)',...
    sprintf('%1.2f (%1.2f-%1.2f)',median_range_ppv(1),median_range_ppv(2),...
    median_range_ppv(3))};
sleep_str = {'Percentage of times classified as asleep: median (range)',...
    sprintf('%1.1f%% (%1.1f-%1.1f)',median_range_sleep(1),median_range_sleep(2),...
    median_range_sleep(3))};



all = [total_str;...
    female_str;...
    age_onset_str;...
    age_implant_str;...
    lat_str;...
    left_str;...
    right_str;...
    bilateral_str;...
    loc_str;...
    temporal_str;...
    other_str;...
    nelecs_str;...
    implant_str;...
    n_gs_str;...
    n_stereo_str;...
    duration_str;...
    ppv_str;...
    rate_str;...
    sleep_str];

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


function any_grid = decide_if_any_grids_or_strips(labels)

any_grid = 0;

for i = 1:length(labels)
    curr = labels{i};
    if contains(curr,'G')
        B = regexp(curr,'\d*','Match');
        if isempty(B), continue; end
        B = str2num(B{1});
        if B > 12 % if it's *G13 or higher then it's a grid
            any_grid = 1;
            break
        end
    end
end

end