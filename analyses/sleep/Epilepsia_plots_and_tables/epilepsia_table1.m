function T1 = epilepsia_table1

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/epilepsia/'];
int_folder = [results_folder,'analysis/intermediate_epilepsia_revision/'];
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
%elec_locs = out.circ_out.all_elec_locs;
%elec_lats = out.circ_out.all_elec_lats;
%[loc_lat_count,patient_locs] = count_locs_and_lats(locs,lats);

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');

%% Load pt file
pt = load([out_folder1,'pt.mat']);
pt = pt.pt;

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
mri_lesional = cell(npts,1);
concordant_loc = cell(npts,1);
concordant_lat = cell(npts,1);

n_one_minute_segments = nan(npts,1);
n_spikes = nan(npts,1);
n_spikes_per_elec = nan(npts,1);
n_seizures = nan(npts,1);

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
    sz_times = summ.sz_times;
    spikes = summ.spikes;
    
    nelecs(p) = length(labels);
    %any_grids(p) = decide_if_any_grids_or_strips(labels);

    %% Basic numbers
    n_one_minute_segments(p) = size(spikes,2);
    n_spikes(p) = nansum(spikes,'all');
    n_spikes_per_elec(p) = nanmean(nansum(spikes,2));
    n_seizures(p) = size(sz_times,1);
    
    %% Spike rate
    rate(p) = nanmean(spikes,'all');
    
    
    %% Duration
    duration(p) = summ.times(end)/3600/24;
    
    %% SOZ localization
    [curr_loc,curr_lat] = seizure_localization_parser(summ.soz.loc,summ.soz.lat);
    lat{p} = curr_lat;
    loc{p} = curr_loc;
    
     %% Spike PPV
    row = strcmp(val_T.name,names{p});
    assert(sum(row) == 1)
    ppv(p) = val_T.PPV_car_(row);
    
    %% Sleep wake PPV
    row = strcmp(sw_val_T.name,names{p});
    assert(sum(row) == 1)
    ws_ppv(p,:) = [sw_val_T.x_Correct_outOf50_Wake(row)/50, sw_val_T.x_Correct_outOf50_Sleep(row)/50];
    
    %% Get preimplant data
    % Find the index in pt structure
    found_it = 0;
    for ip = 1:length(pt)
        if strcmp(summ.name,pt(ip).name)
            found_it = 1;
            mri_lesional{p} = pt(ip).clinical.pre_implant.MRI_lesional;
            concordant_loc{p} = pt(ip).clinical.pre_implant.concordant_loc;
            concordant_lat{p} = pt(ip).clinical.pre_implant.concordant_lat;
             
            break
        end
        
    end
    if found_it == 0, error('why'); end
    
    
end


mri_lesional = cellfun(@(x) clean_preimplant_designations(x),mri_lesional);
concordant_loc = cellfun(@(x) clean_preimplant_designations(x),concordant_loc);
concordant_lat = cellfun(@(x) clean_preimplant_designations(x),concordant_lat);

%{
T1 = table(names,sex,age_onset,age_implant,nelecs,...
    duration,rate,loc,lat);
%}
T1 = table(names,lat,stereo);

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
median_range_nseizures = [nanmedian(n_seizures),min(n_seizures),max(n_seizures)];
n_left = sum(cellfun(@(x) strcmp(x,'left'),lat));
n_right = sum(cellfun(@(x) strcmp(x,'right'),lat));
n_bilateral = sum(cellfun(@(x) strcmp(x,'bilateral'),lat));
n_temporal = sum(cellfun(@(x) strcmp(x,'temporal'),loc));
n_other = sum(cellfun(@(x) strcmp(x,'other'),loc));
n_lesional = sum(mri_lesional==1);
n_nonlesional = sum(mri_lesional==0);
n_concordant_loc = sum(concordant_loc == 1);
n_nonconcordant_loc = sum(concordant_loc == 0);
n_concordant_lat = sum(concordant_lat == 1);
n_nonconcordant_lat = sum(concordant_lat == 0);

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
n_sz_str = {'Number of seizures: median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_nseizures(1),median_range_nseizures(2),median_range_nseizures(3))};
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
mri_str = {'MRI lesional',''};
mri_lesional_str = {'Yes: N (%)',sprintf('%d (%1.1f%%)',n_lesional,n_lesional/npts*100)};
mri_nonlesional_str = {'No: N (%)',sprintf('%d (%1.1f%%)',n_nonlesional,n_nonlesional/npts*100)};
loc_preimplant_str = {'Pre-implant localization concordant',''};
loc_concordant_str = {'Yes: N (%)',sprintf('%d (%1.1f%%)',n_concordant_loc,n_concordant_loc/npts*100)};
loc_nonconcordant_str = {'No: N (%)',sprintf('%d (%1.1f%%)',n_nonconcordant_loc,n_nonconcordant_loc/npts*100)};
lat_preimplant_str = {'Pre-implant lateralization concordant',''};
lat_concordant_str = {'Yes: N (%)',sprintf('%d (%1.1f%%)',n_concordant_lat,n_concordant_lat/npts*100)};
lat_nonconcordant_str = {'No: N (%)',sprintf('%d (%1.1f%%)',n_nonconcordant_lat,n_nonconcordant_lat/npts*100)};


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
    n_sz_str;...
    sleep_str;...
    mri_str;...
    mri_lesional_str;...
    mri_nonlesional_str;...
    loc_preimplant_str;...
    loc_concordant_str;...
    loc_nonconcordant_str;...
    lat_preimplant_str;...
    lat_concordant_str;...
    lat_nonconcordant_str];

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

%% Get some summary stats for text
fprintf(fid,'<p><b>Clinical information, summary of intracranial recording, and sleep staging</b><br>');

fprintf(fid,['<p>18 out of 119 patients were excluded because the positive '...
    'predictive value of their spike detections was less than 70%%. The remaining '...
    '101 patients had predominantly stereo-EEG implantation (82.2%%) (Table 1). '...
    'Across the 101 patients, we studied a total of %d minutes (%1.1f days downsampled from %1.1f days of total recording) of EEG data, '...
    'containing a total of %d spikes across all electrodes, '...
    'and %d total seizures. '...
    'The AUC comparing predicted sleep vs. wake classifications against manual '...
    'sleep vs. wake labels was 0.90 (Fig. 1D). Examining a sample of 50 random '...
    'spike detections from the wake- and sleep-classified periods from '...
    'each patient demonstrated that the spike detector''s positive predictive '...
    'value was similar in wake (median 92%%) and sleep (median 94%%) time periods '...
    '(Wilcoxon sign rank test: <i>T<sup>+</sup></i> = 1439.5, <i>p</i> = 0.12).'],...
    nansum(n_one_minute_segments),nansum(n_one_minute_segments)/60/24,nansum(n_one_minute_segments)/6/24,...
    nansum(n_spikes),nansum(n_seizures));

fclose(fid);

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