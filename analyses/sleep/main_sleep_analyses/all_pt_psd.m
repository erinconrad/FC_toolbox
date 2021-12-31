function out = all_pt_psd

%% Summary
%{
This function calculates the FFT of spike rates to see how circadian spikes
are. It also looks at spike rates as a function of time of day.

It also does some vestigial code to compare the relative circadian
power by localization and spike location.
%}

%% Parameters
main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
main_lats = {'Left','Right'};
main{1} = main_locs;
main{2} = main_lats;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
int_folder = [results_folder,'analysis/intermediate/'];
out_folder = [results_folder,'analysis/sleep/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% Get the longest run (will pad the others with zeros)
longest_run = 0;
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    run_length = length(summ.times);
    if run_length > longest_run
        longest_run = run_length;
    end
end

all_P = cell(2,1);
circ_P = cell(2,1);
for i = 1:length(all_P)
    all_P{i} = nan(length(main{i}),npts,ceil(longest_run/2));
    circ_P{i} = nan(length(main{i}),npts);
end

% Initialize psd
all_psd = nan(npts,ceil(longest_run/2));
fs = 0.0017;%1/summ(1).block_dur;
all_freqs = nan(npts,ceil(longest_run/2));
skip_pts = [];
all_locs = cell(npts,1);
all_lats = cell(npts,1);
all_circ_P = nan(npts,1);
names = cell(npts,1);
sex = cell(npts,1);
age_onset = nan(npts,1);
age_implant = nan(npts,1);
duration = nan(npts,1);

%% for mod midnight, get file size
[~,n_tod_bins,tod_edges] = bin_mod_midnight_times(zeros(5000,1),[]);
all_tod_rate = nan(npts,n_tod_bins); %w, s
%eleven_to_five = nan(npts,2);


%% Loop over patients and get psd per pt
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get patient info
    labels = summ.labels;
    spikes = summ.spikes;
    times = summ.times;
    loc = summ.ana_loc;
    lat = summ.ana_lat;
    run_length = length(times);
    soz_loc = summ.soz.loc;
    soz_lat = summ.soz.lat;
    name = summ.name;
    sex{p} = summ.clinical.sex;
    age_implant(p) = summ.clinical.age_implant;
    age_onset(p) = summ.clinical.age_onset;
    duration(p) = age_implant(p)-age_onset(p);
    mod_midnight = summ.mod_midnight;
    
    %% Find and remove non-intracranial electrodes (ekg and scalp)
    ekg = find_non_intracranial(labels);
    spikes = spikes(~ekg,:);
    loc = loc(~ekg,:);
    lat = lat(~ekg);
    
    %% Get spike rate by time of day
     % Bin the mod midnights
    [mod_midnight] = bin_mod_midnight_times(mod_midnight,tod_edges); % says which of the 144 bins the current time is in
    tod_rate = nan(n_tod_bins,1);
    
    % Loops over bins
    for t = 1:n_tod_bins
        
        % find the times in that time of day bin
        curr_bins = mod_midnight == t; % which runs match that time of day
        tod_rate(t,:) = nansum(spikes(:,curr_bins),'all'); % sum up all spikes across all electrodes and all times in that bin
          
    end
    all_tod_rate(p,:) = tod_rate;
    %eleven_pm_to_five_am = summ.mod_midnight < 5*3600 | summ.mod_midnight > 23*3600;
    %eleven_to_five(p,:) = [nanmean(spikes(:,eleven_pm_to_five_am),'all'),...
    %nanmean(spikes(:,~eleven_pm_to_five_am),'all')];
    names{p} = name;
    
   
    
    % parse SOZ localization
    %[soz_loc,soz_lat] = seizure_localization_parser(soz_loc,soz_lat);
    all_locs{p} = soz_loc;
    all_lats{p} = soz_lat;
    
    % pad spikes
    spikes = [spikes,zeros(size(spikes,1),longest_run-run_length)];
    
    % average across electrodes
    avg_spikes = nanmean(spikes,1);
    
    % get psd
    [P,freqs] = power_by_freq(avg_spikes,fs);
    
    all_psd(p,:) = P;
    all_freqs(p,:) = freqs;
    
    % Get total circular power  
    all_circ_P(p) = get_circ_power(avg_spikes,fs);
    
    
    
    %% Get spectral power for each group for locs and lats
    % Skip if all empty
    if sum(cellfun(@(x) isempty(x),loc)) == length(loc) 
        fprintf('\nskipping pt %d\n',p);
        skip_pts = [skip_pts;p];
        continue
    end
    
    % Loop over loc vs lat
    for g = 1:2
        if g == 1
            group = loc;
        elseif g == 2
            group = lat;
        end
        
        % Get the rates corresponding to the subgroups
        % (can probably do this without a for loop)
       
        for sg = 1:length(main{g})
            ic = ismember(group,main{g}(sg));
            rate_subgroup = nanmean(spikes(ic,:),1);
                    
            % Get the spectral power around 24 hour peak
            all_P{g}(sg,p,:) = power_by_freq(rate_subgroup,fs);
            circ_P{g}(sg,p) = get_circ_power(rate_subgroup,fs);

        end
        

    end
end

%% confirm freqs same across patients
sum_diff_freqs = sum(abs(sum(abs(diff(all_freqs,1,1)))));
assert(sum_diff_freqs < 1e-3);

%% Stuff for PSD
periods = 1./freqs/3600;
low_period = periods <= 100;
periods = periods(low_period);

% Restrict to less than 100
% Get stats
all_psd = all_psd(:,low_period);
all_psd = all_psd./sum(all_psd,2);
median_psd = median(all_psd,1);
iqr_psd = [prctile(all_psd,25,1);prctile(all_psd,75,1)];

%% Prep output
out.all_psd = all_psd;
out.low_period = low_period;
out.median_psd = median_psd;
out.iqr_psd = iqr_psd;
out.periods = periods;
out.circ_P = circ_P;
out.main_locs = main_locs;
out.all_circ_P = all_circ_P;
out.all_locs = all_locs;
out.all_lats = all_lats;
out.skip_pts = skip_pts;
out.names = names;
out.sex = sex;
out.age_onset = age_onset;
out.age_implant = age_implant;
out.duration = duration;
out.all_tod_rate = all_tod_rate;
out.tod_edges = tod_edges;
%out.eleven_to_five = eleven_to_five;

%{
%% Initialize figure
figure
set(gcf,'position',[100 100 1400 500])
tiledlayout(1,4,'tilespacing','tight','padding','tight')

%% Overall spike rate
nexttile
% Plot
mp = shaded_error_bars(periods,median_psd,iqr_psd,[0 0 0]);
xlim([0 100])
xlabel('Period (hours)')
ylabel({'Spike rate', 'normalized power spectrum'});
%legend(mp,'Overall spike rate','fontsize',15,'location','northwest')
set(gca,'fontsize',15)
title('Spike rate power spectral density')

%
%% Do the localizations
nexttile
g = 1;
loc_P = all_P{g};
plotsg = nan(size(loc_P,1),1);
line_type = {'-','--','-.',':'};
for sg = 1:size(loc_P,1)

    % get stats
    curr_psd = squeeze(loc_P(sg,:,:));
    curr_psd = curr_psd(:,low_period);
    curr_psd = curr_psd./sum(curr_psd,2);
    median_psd = nanmedian(curr_psd,1);
    
    %iqr_psd = [prctile(curr_psd,25,1);prctile(curr_psd,75,1)];

    % Plot
    %shaded_error_bars(periods,median_psd,iqr_psd,[]);
    plotsg(sg) = plot(periods,median_psd,line_type{sg},'linewidth',2);
    hold on
 
end
xlim([0 100])
xlabel('Period (hours)')
ylabel({'Spike rate', 'normalized power spectrum'});
set(gca,'fontsize',15)
legend(main_locs)
title({'Spike rate power spectral density','by spike location'});


%% Compare cyclical power across electrode localizations
nexttile
plot_paired_data(circ_P{g},main_locs,'Relative circadian power','paired')
title({'Relative circadian power','by spike location'})

%% Compare cyclical power between patients with different soz localizations
mt = cellfun(@(x) strcmp(x,'temporal'), all_locs);
circ_mt = all_circ_P(mt);
circ_other = all_circ_P(~mt);

% pad both with nans
circ_mt = [circ_mt;nan(length(mt)-length(circ_mt),1)];
circ_other = [circ_other;nan(length(mt)-length(circ_other),1)];
nexttile
plot_paired_data(([circ_mt,circ_other])',{'Temporal','Other'},'Relative circadian power','unpaired')
title({'Relative circadian power','by SOZ localization'});



print([out_folder,'circadian'],'-dpng')
close(gcf)
%}

end

