function sleep_table2

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/neurology/'];
%int_folder = [results_folder,'analysis/intermediate/'];
int_folder = [results_folder,'analysis/backup_intermediate_Feb26_good_spikes/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
out_folder1 = [scripts_folder,'analyses/sleep/data/'];

%% Load validation file
val_T = readtable(['Manual validation.xlsx']);
sw_val_T = readtable(['Manual validation.xlsx'],'Sheet','Validation Sw');

%% Load out file
out = load([out_folder1,'out.mat']);
out = out.out;
n_sleep_wake = out.bin_out.n_sleep_wake;

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

names = cell(npts,1);
nelecs = nan(npts,1);
rate = nan(npts,1);
duration = nan(npts,1);
stereo = nan(npts,1);
ppv = nan(npts,1);
perc_asleep = n_sleep_wake(:,1)./sum(n_sleep_wake,2)*100;
ws_ppv = nan(npts,2);

%% Loop over files
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    name = summ.name;
    
    
    names{p} = summ.name;
    clinical = summ.clinical;
    stereo(p) = clinical.stereo;
    
    %% Electrodes
    labels = summ.labels;
    ekg = find_non_intracranial(labels);
    labels = labels(~ekg);
    
    nelecs(p) = length(labels);
    
    %% Spike rate
    spikes = summ.spikes;
    rate(p) = nanmean(spikes,'all');
    
    %% Duration
    duration(p) = summ.times(end)/3600/24;
    
    %% Spike PPV
    row = strcmp(val_T.Var2,name);
    assert(sum(row) == 1)
    ppv(p) = val_T.Var11(row);
    
    %% Sleep wake PPV
    row = strcmp(sw_val_T.name,name);
    assert(sum(row) == 1)
    ws_ppv(p,:) = [sw_val_T.x_Correct_outOf50_Wake(row)/50, sw_val_T.x_Correct_outOf50_Sleep(row)/50];
    
end

%% Turn into summary stats

% Get data
median_range_nelecs = [nanmedian(nelecs),min(nelecs),max(nelecs)];
median_range_duration = [nanmedian(duration),min(duration),max(duration)];
median_range_rate = [nanmedian(rate),min(rate),max(rate)];
median_range_ppv = [nanmedian(ppv),min(ppv),max(ppv)];
median_range_sleep = [nanmedian(perc_asleep),min(perc_asleep),max(perc_asleep)];
median_range_wake_ppv = [nanmedian(ws_ppv(:,1)),min(ws_ppv(:,1)),max(ws_ppv(:,1))];
median_range_sleep_ppv = [nanmedian(ws_ppv(:,2)),min(ws_ppv(:,2)),max(ws_ppv(:,2))];


% Turn to table
nelecs_str = {'Number of electrodes: median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_nelecs(1),median_range_nelecs(2),median_range_nelecs(3))};
implant_str = {'Implant type',''};
n_gs_str = {'Grids/strips/depths: N (%)',sprintf('%d (%1.1f%%)',sum(stereo==0),sum(stereo==0)/length(stereo)*100)};
n_stereo_str = {'Stereo-EEG: N (%)',sprintf('%d (%1.1f%%)',sum(stereo==1),sum(stereo==1)/length(stereo)*100)};
duration_str = {'Intracranial recording duration in days: median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_duration(1),median_range_duration(2),median_range_duration(3))};
rate_str = {'Spike rate (spikes/elecs/min): median (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    median_range_rate(1),median_range_rate(2),median_range_rate(3))};
ppv_str = {'Spike detector positive predictive value: median (range)',...
    sprintf('%1.2f (%1.2f-%1.2f)',median_range_ppv(1),median_range_ppv(2),...
    median_range_ppv(3))};
sleep_str = {'Percentage of times classified as asleep: median (range)',...
    sprintf('%1.1f%% (%1.1f-%1.1f)',median_range_sleep(1),median_range_sleep(2),...
    median_range_sleep(3))};

all = [nelecs_str;...
    implant_str;...
    n_gs_str;...
    n_stereo_str;...
    duration_str;...
    ppv_str;...
    rate_str;...
    sleep_str];

T2 = cell2table(all);
writetable(T2,[out_folder,'Table2.csv']);

%% Results sentence comparing ws ppv
[pval,~,stats] = signrank(ws_ppv(:,1),ws_ppv(:,2));
Tpos =stats.signedrank;
fprintf(['\nExamining a sample of 50 random spike detections from '...
    'the wake- and sleep-classified periods from each patient demonstrated '...
    'that the spike detector PPV was similar in wake (median PPV %1.1f%%) '...
    'and sleep (median PPV %1.1f%%) time periods (T+ = %1.1f, '...
    'p = %1.2f).\n'],nanmedian(ws_ppv(:,1))*100,nanmedian(ws_ppv(:,2))*100,Tpos,pval);

end