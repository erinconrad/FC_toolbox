function ccep_pc_dist(pc,ccep)

% change wrap method

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

%% Get pc info
out = net_over_time(pc);
out = reconcile_files(out);

%% Get all base labels
% pc
base_pc_labels = out.all_labels;

% ccep
ccep_labels = ccep.chLabels;


%% Find intersecting labels and corresponding indices
% we need to remember to restrict all future things based on these indices
[indices,intersecting_set] = find_intersecting_idx({base_pc_labels,ccep_labels});
pc_idx = indices{1};
ccep_idx = indices{2};

%% Get avg pc network
pc_net = out.montage(m).net;

% Take average over all times for net
pc_net = nanmean(pc_net,2);

% unwrap
pc_net= wrap_or_unwrap_adjacency_fc_toolbox(pc_net);
pc_net(logical(eye(size(pc_net)))) = nan;
nan_chs = sum(~isnan(pc_net),1) == 0;

% restrict indices
pc_labels = out.montage(m).labels(pc_idx);
pc_net = pc_net(pc_idx,pc_idx);

% remove nans
%{
pc_net(nan_chs,:) = [];
pc_net(:,nan_chs) = [];
pc_labels(nan_chs) = [];
%}

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

if m == 1
    ccep_labels = ccep_labels_bipolar;
elseif m == 2
    ccep_labels = ccep_labels_car;
end

ccep_labels = ccep_labels(ccep_idx);
ccep_net = ccep_net(ccep_idx,ccep_idx);
%{
stim_labels = ccep_labels(stim_chs);
response_labels = ccep_labels(response_chs);

% Restrict matrix
ccep_net = ccep_net(response_chs,stim_chs);

% calculate indegree and outdegree
outdegree = nansum(ccep_net,1); % returns 1 x nchs column vector, 1 for each stim (but need to reduce chs!)
indegree = nansum(ccep_net,2); % returns nchs x 1 row vector
outdegree(sum(~isnan(ccep_net),1)==0) = nan;
indegree(sum(~isnan(ccep_net),2)==0) = nan;

% Make non-stim or non-response chs nans
outdegree(~stim_chs) = nan;
outdegree = outdegree';
indegree(~response_chs) = nan;
%}

%% Get distance stuff
% Reconcile locs and anatomy with ieeg labels
locs = pt(p).elecs(1).locs;
loc_labels = decompose_labels(pt(p).elecs(1).elec_names);
locs = reconcile_locs_ieeg(base_pc_labels,loc_labels,locs,[]);

% locs
[~,~,bipolar_labels,chs_in_bipolar,~,mid_locs,~] =...
    bipolar_montage(nan(1,length(base_pc_labels)),base_pc_labels,1:length(base_pc_labels),locs,[]);
% Made it to here!!!!!!!!!!!

[~,car_labels] = car_montage(nan(1,length(base_pc_labels)),1:length(base_pc_labels),base_pc_labels);

if m ==1
    labels = bipolar_labels;
    out_locs = mid_locs;
else
    out_locs = locs;
    labels = car_labels;
end

% restrict idx
out_locs = out_locs(pc_idx,:);
labels = labels(pc_idx);


% Make inv dist matrix
A = inverse_dist(out_locs);
A = wrap_or_unwrap_adjacency_fc_toolbox(A);
%{
A(nan_chs,:) = [];
A(:,nan_chs) = [];
dist_labels(nan_chs) = [];
%}

%% Find common labels and restrict
%{
% common labels
stim_overlap_labels = intersect(dist_labels,intersect(stim_labels,pc_labels));
response_overlap_labels = intersect(dist_labels,intersect(response_labels,pc_labels));

% restrict stim
stim_stim_idx = ismember(stim_labels,stim_overlap_labels);
ccep_net(:,~stim_stim_idx) = [];
outdegree(~stim_stim_idx) = [];
stim_labels(~stim_stim_idx) = [];
%}

%% Reconcile labels
indices = reconcile_labels({ccep_labels,pc_labels,labels});

% Restrict matrices to matching labels
ccep_labels = ccep_labels(indices{1});
ccep_net = ccep_net(indices{1},indices{1});

pc_labels = pc_labels(indices{2});
pc_net = pc_net(indices{2},indices{2});

labels = labels(indices{3});
A = A(indices{3},indices{3});


%% Make sure labels match
%
if ~isequal(pc_labels,ccep_labels) || ~isequal(pc_labels,labels)
    error('labels do not match')
end
%}

%% Do correlations
% Correlations across columns
outdegree = nansum(ccep_net,1);
ns_pc_cols = nansum(pc_net,1);
ns_dist_cols = nansum(A,1);
out_pc_r = corr(outdegree',ns_pc_cols','rows','pairwise');
out_dist_r = corr(outdegree',ns_dist_cols','rows','pairwise');
pc_dist_r_cols = corr(ns_dist_cols',ns_pc_cols','rows','pairwise');

% Correlations across rows
indegree = nansum(ccep_net,2);
ns_pc_rows = nansum(pc_net,2);
ns_dist_rows = nansum(A,2);
in_pc_r = corr(indegree,ns_pc_rows,'rows','pairwise');
in_dist_r = corr(indegree,ns_dist_rows,'rows','pairwise');
pc_dist_r_rows = corr(ns_dist_rows,ns_pc_rows,'rows','pairwise');

if abs(pc_dist_r_rows - pc_dist_r_cols) > 1e-4, error('what'); end
fprintf('\nDistance-PC r = %1.2f\n',pc_dist_r_rows);
fprintf('\nOutdegree-distance r = %1.2f\n',out_dist_r);
fprintf('\nOutdegree-PC r = %1.2f\n',out_pc_r);
fprintf('\nIndegree-distance r = %1.2f\n',in_dist_r);
fprintf('\nIndegree-PC r = %1.2f\n',in_pc_r);


%% Plot
figure
set(gcf,'position',[10 10 1300 900])
tiledlayout(2,3,'tilespacing','tight','padding','tight')

nexttile
turn_nans_gray(A)
xticks(1:length(labels))
yticks(1:length(labels))
xticklabels(labels)
yticklabels(labels)
title('Inverse distance')

nexttile
turn_nans_gray(pc_net)
xticks(1:length(pc_labels))
yticks(1:length(pc_labels))
xticklabels(pc_labels)
yticklabels(pc_labels)
title('Pearson correlation')

nexttile
turn_nans_gray(ccep_net)
xticks(1:length(ccep_labels))
yticks(1:length(ccep_labels))
xticklabels(ccep_labels)
yticklabels(ccep_labels)
xlabel('Stim')
ylabel('Response')
title('CCEP')

nexttile
plot(ns_dist_rows,ns_pc_rows,'o')
yl = ylim;
xl = xlim;
text(xl(1),yl(2),sprintf('r = %1.2f\n',pc_dist_r_rows),'verticalalignment','top',...
    'fontsize',15)
xlabel('Inverse distance node strength')
ylabel('PC node strength')

nexttile
plot(ns_pc_rows,indegree,'o')
yl = ylim;
xl = xlim;
text(xl(1),yl(2),sprintf('r = %1.2f\n',in_dist_r),'verticalalignment','top',...
    'fontsize',15)
xlabel('PC node strength')
ylabel('Indegree')

nexttile
plot(ns_pc_rows,outdegree,'o')
yl = ylim;
xl = xlim;
text(xl(1),yl(2),sprintf('r = %1.2f\n',out_dist_r),'verticalalignment','top',...
    'fontsize',15)
xlabel('PC node strength')
ylabel('Outdegree')



end