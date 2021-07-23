function ccep_pc_dist(pc,ccep)

%% Parameters
m = 1;

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Pt struct
data_folder = [locations.main_folder,'data/'];
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% find corresponding p in pt
name = pc.name;
for p = 1:length(pt)
    if strcmp(name,pt(p).name)
        break
    end
end

%% Get avg pc network
out = net_over_time(pc);
out = reconcile_files(out);
pc_labels = out.montage(m).labels;
pc_net = out.montage(m).net;

% Take average over all times for net
pc_net = nanmean(pc_net,2);

% unwrap
pc_net= wrap_or_unwrap_adjacency_fc_toolbox(pc_net);
pc_net(logical(eye(size(pc_net)))) = nan;
nan_chs = sum(~isnan(pc_net),1) == 0;

% remove nans for plotting
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

%% Get distance stuff
% Reconcile locs and anatomy with ieeg labels
loc_labels = pt(p).elecs.elec_names;
locs = pt(p).elecs.locs;
anatomy = pt(p).elecs.anatomy;
clean_loc_labels = decompose_labels(loc_labels);
clean_labels = out.all_labels;
[locs,anatomy] = reconcile_locs_ieeg(clean_labels,clean_loc_labels,locs,anatomy);

% locs
[~,~,bipolar_labels,chs_in_bipolar,~,mid_locs,~] =...
    bipolar_montage(nan(1,length(clean_labels)),clean_labels,1:length(clean_labels),locs,anatomy);

[~,car_labels] = car_montage(nan(1,length(clean_labels)),1:length(clean_labels),clean_labels);

if m ==1
    labels = bipolar_labels;
    out_locs = mid_locs;
else
    out_locs = locs;
    labels = car_labels;
    
end


% Make inv dist matrix
A = inverse_dist(out_locs);
A = wrap_or_unwrap_adjacency_fc_toolbox(A);
A(nan_chs,:) = [];
A(:,nan_chs) = [];
labels(nan_chs) = [];

%% Make sure labels match
%
if ~isequal(pc_labels,ccep_labels) || ~isequal(pc_labels,labels)
    error('labels do not match')
end
%}

%% Plot
figure
set(gcf,'position',[10 10 1300 500])
tiledlayout(1,3,'tilespacing','tight','padding','tight')

nexttile
turn_nans_gray(pc_net)
xticks(1:length(pc_plot_labels))
yticks(1:length(pc_plot_labels))
xticklabels(pc_plot_labels)
yticklabels(pc_plot_labels)

nexttile
turn_nans_gray(A)
xticks(1:length(labels))
yticks(1:length(labels))
xticklabels(labels)
yticklabels(labels)

nexttile
turn_nans_gray(ccep_net)
xticks(1:length(ccep_stim_labels))
yticks(1:length(ccep_response_labels))
xticklabels(ccep_stim_labels)
yticklabels(ccep_response_labels)


end