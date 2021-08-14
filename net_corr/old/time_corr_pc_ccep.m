function time_corr_pc_ccep(pc,ccep)

%% Parameters
f = 1;
m = 1;

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Get ns
pc_ns_struct = ns_over_time(pc);
pc_ns = pc_ns_struct.file(f).montage(m).data;
pc_labels = pc.file(f).run(1).data.montage(m).labels;
pc_net = pc_ns_struct.file(f).montage(m).net;
pc_net(logical(eye(size(pc_net)))) = nan;
nan_chs = sum(~isnan(pc_net),1) == 0;
pc_net(nan_chs,:) = [];
pc_net(:,nan_chs) = [];
pc_plot_labels = pc_labels;
pc_plot_labels(nan_chs) = [];

%% Get ccep stuff
% network and labels
ccep_net = ccep.A;
ccep_labels_bipolar = ccep.bipolar_labels;
ccep_labels_car = ccep.chLabels;

% put labels into same format as pc labels
empty_labels = cellfun(@isempty,ccep_labels_bipolar);
ccep_labels_bipolar(empty_labels) = {'-'};
ccep_labels_car = cellfun(@(x) [x,'-CAR'],ccep_labels_car,'UniformOutput',false);

% Identify stim and response chs
stim_chs = ccep.ch_info.stim_chs;
response_chs = ccep.ch_info.response_chs;

% Restrict matrix
ccep_net = ccep_net(response_chs,stim_chs);

if m == 1
    ccep_labels = ccep_labels_bipolar;
elseif m == 2
    ccep_labels = ccep_labels_car;
end

ccep_stim_labels = ccep_labels(stim_chs);
ccep_response_labels = ccep_labels(response_chs);

% calculate indegree and outdegree
outdegree = nansum(ccep_net,1); % returns 1 x nchs column vector, 1 for each stim (but need to reduce chs!)
indegree = nansum(ccep_net,2); % returns nchs x 1 row vector
outdegree(sum(~isnan(ccep_net),1)==0) = nan;
indegree(sum(~isnan(ccep_net),2)==0) = nan;

% Make non-stim or non-response chs nans
outdegree(~stim_chs) = nan;
outdegree = outdegree';
indegree(~response_chs) = nan;



%% Make sure labels match
if ~isequal(pc_labels,ccep_labels) 
    error('labels do not match')
end

%% Get pairwise correlations for each time
[ns_out_corr,out_p] = corr(pc_ns,outdegree,'rows','pairwise');
[ns_in_corr,in_p] = corr(pc_ns,indegree,'rows','pairwise');



%% Plots
figure
set(gcf,'position',[10 10 1200 1000])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

nexttile
turn_nans_gray(pc_net)
xticks(1:length(pc_plot_labels))
yticks(1:length(pc_plot_labels))
xticklabels(pc_plot_labels)
yticklabels(pc_plot_labels)

nexttile
turn_nans_gray(ccep_net)
xticks(1:length(ccep_stim_labels))
xticklabels(ccep_stim_labels)
yticks(1:length(ccep_response_labels))
yticklabels(ccep_response_labels)

nexttile([1,2])
plot(ns_out_corr)
hold on
plot(ns_in_corr)
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('NS-out r = %1.2f\nNS-in r = %1.2f',nanmean(ns_out_corr),nanmean(ns_in_corr)),...
    'verticalalignment','top')
legend({'NS-outdegree','NS-indegree'})

end