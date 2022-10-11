function out = sleep_seeg_ad

% Cool! It looks like the normalized AD ratio agrees pretty well with the
% Sleep SEEG designations if I simplify them to wake/sleep. Interestingly,
% REM sleep also has a very low AD ratio (though not as low as N3).


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_out_dir = [results_folder,'edf_out/'];
sleep_manual_dir = [results_folder,'analysis/sleep/erin_designations/'];
int_folder = [results_folder,'analysis/backup_intermediate_Feb26_good_spikes/'];
out_file = [results_folder,'analysis/sleep/sleep_seeg/'];
if ~exist(out_file,'dir')
    mkdir(out_file)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

all_out = {};

%% Loop over patients
for p = 1:npts
    
    pt_name = strrep(listing(p).name,'.mat','');
    if ~exist([edf_out_dir,pt_name,'/sleep_stage.mat'])
        continue
    end

    fprintf('\nDoing patient %s\n',pt_name);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;

    %% Load sleepseeg designation file
    sout = load([edf_out_dir,pt_name,'/sleep_stage.mat']);
    sout = sout.sout;

    % Get times of sleep transitions
    st_dates= sout.Summary(2:end,2);
    st_times = sout.Summary(2:end,3);
    stage = sout.Summary(2:end,4);
    
    % convert to seconds in ieeg
    seeg_secs = convert_transitions_to_ieeg_secs(st_dates,st_times);

    %% Get main things
    labels = summ.labels;
    ad = summ.ad;
    file_index = summ.file_index;
    times = summ.times;

    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    ad_norm = (ad-nanmedian(ad))./iqr(ad);

    % get the sleep stages for appropriate times
    [is_seeg_time,seeg_stage] = match_seeg_times(times,file_index,seeg_secs,stage);

    % fill out all out
    all_out = [all_out;...
        repmat({pt_name},sum(is_seeg_time),1),seeg_stage(is_seeg_time),num2cell(ad_norm(is_seeg_time)')];
end


%% Simple labels
ntimes = size(all_out,1);
labels = cell(ntimes,1);
labels(cellfun(@(x) strcmp(x,'W'),all_out(:,2))) = {'Wake'};
labels(cellfun(@(x) ismember(x,{'R','N1','N2','N3'}),all_out(:,2))) = {'Sleep'};

alt_labels = cell(ntimes,1);
alt_labels(cellfun(@(x) strcmp(x,'W'),all_out(:,2))) = {'Wake'};
alt_labels(cellfun(@(x) ismember(x,{'N3'}),all_out(:,2))) = {'Sleep'};
empty_labels = cellfun(@isempty,alt_labels);

%% Scores
scores = cell2mat(all_out(:,3));

[X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'Wake');
[altX,altY,~,altAUC] = perfcurve(alt_labels(~empty_labels),scores(~empty_labels),'Wake');

out.X = X;
out.Y = Y;
out.AUC = AUC;
out.altX = altX;
out.altY = altY;
out.altAUC = altAUC;

%% Show scores for different states
if 0
    figure
    set(gcf,'position',[289 517 1001 280])
    tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
    nexttile
    plot(X,Y,'linewidth',2)
    hold on
    plot(altX,altY,'linewidth',2)
    xlabel('False positive rate')
    ylabel('True positive rate')
    set(gca,'fontsize',15)
    title('ROC classifying SleepSEEG sleep by normalized ADR')
    legend({sprintf('All sleep vs wake AUC %1.2f',AUC),sprintf('N3 vs wake AUC %1.2f',altAUC)},'location','southeast')
    

    nexttile
    boxplot(scores,all_out(:,2))
    set(gca,'fontsize',15)
    ylabel('Normalized alpha-delta ratio')
    xlabel('SleepSEEG classification')
    title('Normalized ADR by SleepSEEG classification')

    
    print(gcf,[out_file,'seeg_vs_ad'],'-dpng')
    close gcf
end

