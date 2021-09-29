function hub_stability

%% Parameters
m = 2; % do not change

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/simple/stabilities/'];
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

% Loop over patients
for l = 1:npts
    j = good_pts(l);
    name = pt(j).name;
    
    %% Load the spike file
    fname = [spikes_folder,name,'_pc.mat'];
    
    if ~exist(fname,'file')
        
        fprintf('\n%s unavailable, skipping\n',name);
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
    times = out.times; times = times/3600/25; % convert times to days
    spikes = out.montage(m).spikes;
    labels = out.montage(m).labels;
    run_dur = diff(pt(j).ieeg.file(1).run_times(1,:)); % 60 s
    spikes = spikes/run_dur*60; % convert spikes to spikes/minute (note this divides by 60 and then multiplies by 60 and so does nothing)
    net = out.montage(m).net;
    
    % Calculate the ns
    net = wrap_or_unwrap_adjacency_fc_toolbox(net);
    ns = squeeze(nansum(net,1));
    
    
    ekg = find_non_intracranial(labels);
    % remove non intracranial
    labels(ekg) = [];
    spikes(ekg,:) = [];
    ns(ekg,:) = []; 

    spikes = spikes/length(labels); % get # spikes/elec
    
    % Clean the labels
    labels = decompose_labels(labels,name);    
    
    % Calculate the average across times
    spike_avg = nanmean(spikes,2);
    ns_avg = nanmean(ns,2);
    
    % calculate the stability over time
    spike_stability = corr(spikes,spike_avg,'Type','Spearman','rows','pairwise');
    ns_stability = corr(ns,ns_avg,'Type','Spearman','rows','pairwise');
    
    % Calculate correlation between spikes and ns
    ntimes = length(times);
    spike_ns_agreement = nan(ntimes,1);
    for it = 1:ntimes
        spike_ns_agreement(it) = corr(spikes(:,it),ns(:,it),'Type','Spearman','rows','pairwise');
    end
    
    % Plot these
    figure
    set(gcf,'position',[10 10 1200 900])
    plot(times,spike_stability,'linewidth',2)
    hold on
    plot(times,ns_stability,'linewidth',2)
    plot(times,spike_ns_agreement,'linewidth',2)
    plot(xlim,[0 0],'k--','linewidth',2)
    xlim([times(1) times(end)])
    xlabel('Day')
    ylim([-1 1])
    title(name)
    legend({'Spike stability','NS stability','Spike-NS agreement'})
    set(gca,'fontsize',15);
    
    % save and close the plot
    print(gcf,[out_folder,name],'-dpng')
    close(gcf)
    
end


end