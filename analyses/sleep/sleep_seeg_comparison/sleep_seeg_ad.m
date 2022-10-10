function sleep_seeg_ad

% Cool! It looks like the normalized AD ratio agrees pretty well with the
% Sleep SEEG designations if I simplify them to wake/sleep. Interestingly,
% REM sleep also has a very low AD ratio (though not as low as N3).

ref_start_time = datetime('01/01/2000 00:00:00','InputFormat','MM/dd/yyyy hh:mm:ss');

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
    
    new_st_dates = st_dates;
    for i = 1:length(new_st_dates) 
        if strcmp(new_st_dates{i},' ') % replace empty with last date
            new_st_dates{i} = new_st_dates{i-1};
        end
    end
    st_dates = new_st_dates;
    
    % Combine dates and times
    dts = cellfun(@(x,y) [x, ' ',y],st_dates,st_times,'uniformoutput',false);

    % convert to seconds into ieeg file
    dt = cellfun(@(x) datetime(x,'InputFormat','dd-MMM-yyyy HH:mm:ss'),dts);
    seeg_secs = seconds(dt-ref_start_time);    
    
    %% Get main things
    labels = summ.labels;
    ad = summ.ad;
    file_index = summ.file_index;
    times = summ.times;

    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    ad_norm = (ad-nanmedian(ad))./iqr(ad);

    %% Loop it
    for t = 1:length(times)
        if file_index(t) > 1
            break
        end

        curr_time = times(t);
        run_times = [curr_time-30,curr_time+30]; % time is middle of 60 second window

         % Figure out if first run time falls between two transitions
        run_start_minus_seeg = run_times(1)-seeg_secs;

        if all(run_start_minus_seeg>0) % this block is too late, stop loop
            break
        end

        if all(run_start_minus_seeg<0) % this block too early, go to next loop
            continue
        end

        % find the first seeg transition that comes BEFORE
        prior_transition = find(run_start_minus_seeg<0);
        prior_transition = prior_transition(1)-1;

        % confirm that run time end comes before the next one
        if run_times(2) > seeg_secs(prior_transition+1)
            continue
        end

        % If made it here, then the run time falls between two state
        % transitions, and so I should be able to compare
        all_out = [all_out;pt_name,stage(prior_transition),ad_norm(t)];
    end
end


%% Simple labels
ntimes = size(all_out,1);
labels = cell(ntimes,1);
labels(cellfun(@(x) strcmp(x,'W'),all_out(:,2))) = {'Wake'};
labels(cellfun(@(x) ismember(x,{'R','N1','N2','N3'}),all_out(:,2))) = {'Sleep'};

%% Scores
scores = cell2mat(all_out(:,3));

[X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'Wake');

%% Show scores for different states
if 1
    figure
    set(gcf,'position',[289 517 1001 280])
    tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
    nexttile
    boxplot(scores,all_out(:,2))
    set(gca,'fontsize',15)
    ylabel('Normalized alpha-delta ratio')
    xlabel('SleepSEEG classification')
    title('Normalized ADR by SleepSEEG classification')

    nexttile
    plot(X,Y,'linewidth',2)
    xlabel('False positive rate')
    ylabel('True positive rate')
    set(gca,'fontsize',15)
    title('ROC classifying SleepSEEG wake vs sleep by normalized ADR')
    legend(sprintf('AUC %1.2f',AUC),'location','southeast')
    print(gcf,[out_file,'seeg_vs_ad'],'-dpng')
    close gcf
end

