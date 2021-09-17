function spikes_cyclic_anatomy

%{
To do:
- add stats
- remove szs
- look at effect of soz
- look at spike leader
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
spikes_folder = [results_folder,'all_out/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

validation_file = [scripts_folder,'spike_detector/Manual validation.xlsx'];

% Pt struct
data_folder = [locations.main_folder,'data/'];
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Get the indices of the patients with good spikes
T = readtable(validation_file);
good_pts = T.Var13;
good_pts = good_pts(~isnan(good_pts));
npts = length(good_pts);


all_rates = cell(2,1);
all_P = cell(2,1);
avg_rate = cell(2,1);
for i = 1:length(all_P)
    avg_rate{i} = nan(length(main{i}),npts);
    all_P{i} = nan(length(main{i}),npts);
end

% Loop over patients
for l = 1:npts
    j = good_pts(l);
    name = pt(j).name;
    
    %% Load the spike file
    fname = [spikes_folder,name,'_pc.mat'];
    
    if ~exist(fname,'file')
        for i = 1:length(all_P)
            avg_rate{i}(:,l) = nan(length(main{i}),1);
            all_P{i}(:,l) = nan(length(main{i}),1);
            fprintf('\n%s unavailable, skipping\n',name);
        end
        continue
    end
    
    %% Get basic info from the patient
    % load the spike file
    pc = load(fname);
    pc = pc.pc;
    
    % Skip the patient if it's incomplete
    if length(pc.file) < length(pt(j).ieeg.file) || ...
            length(pc.file(end).run) < size(pt(j).ieeg.file(end).run_times,1)
        
        fprintf('\n%s incomplete, skipping\n',name);
        continue
    end
        
    
    % reconcile files (deal with changes in electrode names)
    out = net_over_time(pc);
    out = reconcile_files(out);
    
    % Get the spikes and the labels
    times = out.times;
    spikes = out.montage(m).spikes;
    labels = out.montage(m).labels;
    
    % Clean the labels
    clean_labels = decompose_labels(labels,name);    
    
    % Get number of electrode localizations
    ne = length(pt(j).elecs);
    
    %% Find the anatomy corresponding to the spike labelsomy
    % Initialize a new spike_anatomy cell array corresponding to spike
    % labels
    spike_anatomy = cell(length(labels),1);
    spike_locs = nan(length(labels),3);
    already_filled = zeros(length(labels),1);
    
    for e = 1:ne
        % Get loc/anatomy names and labels
        ana_name = pt(j).elecs(e).elec_names;
        locs = pt(j).elecs(e).locs;
        anatomy = pt(j).elecs(e).anatomy;

        % Indices of the loc/anatomy names that match the spike labels
        [lia,locb] = ismember(clean_labels,ana_name);
        % sanity check
        if ~isequal(clean_labels(lia~=0 & already_filled == 0),ana_name(locb(lia~=0 & already_filled == 0))), error('oh no'); end

        
        % Fill up spike anatomy and locs with the anatomy
        if ~strcmp(class(anatomy),'double')
            spike_anatomy(lia~=0 & already_filled == 0) = anatomy(locb(lia~=0 & already_filled == 0));
        end
        spike_locs(lia~=0 & already_filled == 0,:) = locs(locb(lia~=0 & already_filled == 0),:);
        
        % set already filled
        already_filled(lia~=0 & already_filled == 0) = 1;
    end
    
    % WRITE SOMETHING TO CHECK HOW MANY ARE EMPTY 
    
    % Get anatomical groupings
    [loc,lat] = cluster_anatomical_location(spike_anatomy);
    
    
    %% Get rates and spectral power for each group for these groups
    % Loop over loc vs lat
    for g = 1:2
        if g == 1
            group = loc;
        elseif g == 2
            group = lat;
        end
        
        % Get the rates corresponding to the subgroups
        % (can probably do this without a for loop)
        rate_subgroup = nan(length(main{g}),size(spikes,2));
        p_subgroup = nan(length(main{g}),1);
        rate_avg_ind = nan(length(main{g}),1);
        for sg = 1:length(main{g})
            ic = ismember(group,main{g}(sg));
            rate_subgroup(sg,:) = nanmean(spikes(ic,:),1);
                    
            % Get the spectral power around 24 hour peak
            p_subgroup(sg) = get_circ_power(rate_subgroup(sg,:),times);
            
            % Average rate over time
            rate_avg_ind(sg) = nanmean(rate_subgroup(sg,:));
            
        end
        
   
        all_rates{g} = rate_subgroup;
        all_P{g}(:,l) = p_subgroup;
        avg_rate{g}(:,l) = rate_avg_ind;
    end
    
    
    
    %% Plots for single pt
    if 0
        figure
        tiledlayout(2,2)
        
        for g = 1:length(all_rates)
            curr_rates = all_rates{g};
            nexttile
            plot(curr_rates')
            legend(main{g})
            
            nexttile
            plot(all_P{g}(:,l),'o')
            xticks(1:length(all_P{g}(:,l)))
            xticklabels(main{g})
        end
    end
   
    
end


%% Aggregate plots
figure
tiledlayout(2,2)

for g = 1:length(all_rates)
    
    %% Rate
    nexttile
    curr_rate = avg_rate{g}/10; % divide by 10 minutes
    avg_over_pts = nanmean(curr_rate,2);
    std_over_pts = nanstd(curr_rate,[],2);
    ngroups = size(curr_rate,1);
    npts = size(curr_rate,2);
    
    % Do stats
    [p,stats,post_hoc_p,which_groups] = circ_stats(curr_rate);

    errorbar(avg_over_pts,std_over_pts,'o','markersize',10);
    hold on
    xticks(1:length(avg_over_pts))
    xticklabels(main{g})
    xtickangle(30)
    xlim([0 length(avg_over_pts)+1])
    ylabel('Spikes/elec/minute')
    
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
            [heights(size(heights,1)-1,1) heights(size(heights,1)-1,1)],'k');
        text(mean([1 length(avg_over_pts)]),heights(size(heights,1)-1,2),...
            'ns','fontsize',20,'horizontalalignment','center')
    else
        for k = 1:size(pairs_to_plot,1)
            plot([pairs_to_plot(k,1)+0.1 pairs_to_plot(k,2)-0.1],[heights(k,1) heights(k,1)],'k-')
            hold on
            text(mean(pairs_to_plot(k,:)),pairs_to_plot(k,2),...
                get_asterisks(post_hoc_p_to_plot(k),size(which_groups,1)),...
                'fontsize',20,'horizontalalignment','center')
        end
    end
    
    %% power

    nexttile
    curr_power = all_P{g};
    avg_over_pts = nanmean(curr_power,2);
    std_over_pts = nanstd(curr_power,[],2);
    errorbar(avg_over_pts,std_over_pts,'o','markersize',10);
    hold on
    xticks(1:length(avg_over_pts))
    xticklabels(main{g})
    xtickangle(30)
    ylabel('Relative circadian power')
    xlim([0 length(avg_over_pts)+1])
    
    % Do stats
    [p,stats,post_hoc_p,which_groups] = circ_stats(curr_rate);
    
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
            [heights(size(heights,1)-1,1) heights(size(heights,1)-1,1)],'k');
        text(mean([1 length(avg_over_pts)]),heights(size(heights,1)-1,2),...
            'ns','fontsize',10,'horizontalalignment','center')
    else
        for k = 1:size(pairs_to_plot,1)
            plot([pairs_to_plot(k,1)+0.1 pairs_to_plot(k,2)-0.1],[heights(k,1) heights(k,1)],'k-')
            hold on
            text(mean(pairs_to_plot(k,:)),pairs_to_plot(k,2),...
                get_asterisks(post_hoc_p_to_plot(k),size(which_groups,1)),...
                'fontsize',10,'horizontalalignment','center')
        end
    end
end
print([out_folder,'circ_power'],'-dpng');
close(gcf)


end