function all_pt_psd(summ)

%% To do
%{
- add annotations
- add patient-level analyses
%}

%% Parameters
m = 2; % do not change
main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
main_lats = {'Left','Right'};
main{1} = main_locs;
main{2} = main_lats;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Get the longest run (will pad the others with zeros)
longest_run = 0;
npts = length(summ);
for p = 1:npts
    run_length = length(summ(p).times);
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

%% Loop over patients and get psd per pt
for p = 1:npts
    spikes = summ(p).spikes;
    times = summ(p).times;
    loc = summ(p).ana_loc;
    lat = summ(p).ana_lat;
    run_length = length(times);
    
    % pad spikes
    spikes = [spikes,zeros(size(spikes,1),longest_run-run_length)];
    
    % average across electrodes
    avg_spikes = nanmean(spikes,1);
    
    % get psd
    [P,freqs] = power_by_freq(avg_spikes,fs);
    
    all_psd(p,:) = P;
    all_freqs(p,:) = freqs;
    
    %% Get spectral power for each group for locs and lats
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

%% Initialize figure
periods = 1./freqs/3600;
low_period = periods <= 100;
periods = periods(low_period);
figure
set(gcf,'position',[100 100 1100 300])
tiledlayout(1,3,'tilespacing','tight','padding','tight')

%% Overall spike rate
nexttile
% Restrict to less than 100
% Get stats
all_psd = all_psd(:,low_period);
all_psd = all_psd./sum(all_psd,2);
median_psd = median(all_psd,1);
iqr_psd = [prctile(all_psd,25,1);prctile(all_psd,75,1)];

% Plot
mp = shaded_error_bars(periods,median_psd,iqr_psd,[0 0 0]);
xlim([0 100])
xlabel('Period (hours)')
ylabel({'Spike rate', 'normalized power spectrum'});
legend(mp,'Overall spike rate','fontsize',15,'location','northwest')
set(gca,'fontsize',15)

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

%% Compare cyclical power across localizations
nexttile
curr_power = circ_P{g};
median_over_pts = nanmedian(curr_power,2);
iqr_over_pts = [prctile(curr_power,25,2),prctile(curr_power,75,2)];
errorbar(1:length(main_locs),median_over_pts,...
    median_over_pts-iqr_over_pts(:,1),...
    iqr_over_pts(:,2)-median_over_pts,'o','markersize',10,...
    'color','k','linewidth',2);
hold on
xticks(1:length(median_over_pts))
xticklabels(main{g})
xtickangle(30)
ylabel('Relative circadian power')
set(gca,'fontsize',15)
xlim([0 length(median_over_pts)+1])

% Do stats
[p,post_hoc_p,which_groups] = non_para_circ_stats(curr_power);

if p > 0.05
    pairs_to_plot = [];
else
    pairs_to_plot = which_groups(post_hoc_p < 0.05/size(which_groups,1),:);
    post_hoc_p_to_plot = post_hoc_p(post_hoc_p < 0.05/size(which_groups,1));
end
yl=ylim;
heights = get_heights(yl,pairs_to_plot);
ylim([yl(1) heights(end,2)]);

if p > 0.05
    plot([1 length(avg_over_pts)],...
        [heights(size(heights,1)-1,1) heights(size(heights,1)-1,1)],'k',...
        'linewidth',2);
    text(mean([1 length(avg_over_pts)]),heights(size(heights,1)-1,2),...
        'ns','fontsize',15,'horizontalalignment','center')
else
    for k = 1:size(pairs_to_plot,1)
        plot([pairs_to_plot(k,1)+0.1 pairs_to_plot(k,2)-0.1],[heights(k,1) heights(k,1)],'k-',...
            'linewidth',2)
        hold on
        text(mean(pairs_to_plot(k,:)),heights(k,2),...
            get_asterisks(post_hoc_p_to_plot(k),size(which_groups,1)),...
            'fontsize',15,'horizontalalignment','center')
    end
end

print([out_folder,'circadian'],'-dpng')

end

