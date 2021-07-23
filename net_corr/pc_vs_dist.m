function pc_vs_dist(pc)

%% Parameters
f = 1;
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

%% Get ns
pc_ns_struct = ns_over_time(pc);
pc_ns = pc_ns_struct.file(f).montage(m).data;
clean_labels = pc.file(f).run(1).data.clean_labels;
pc_labels = pc.file(f).run(1).data.montage(m).labels;
pc_net = pc_ns_struct.file(f).montage(m).net;
pc_net(logical(eye(size(pc_net)))) = nan;
nan_chs = sum(~isnan(pc_net),1) == 0;
pc_net(nan_chs,:) = [];
pc_net(:,nan_chs) = [];
pc_plot_labels = pc_labels;
pc_plot_labels(nan_chs) = [];

%% Reconcile locs and anatomy with ieeg labels
loc_labels = pt(p).elecs.elec_names;
locs = pt(p).elecs.locs;
anatomy = pt(p).elecs.anatomy;
clean_loc_labels = decompose_labels(loc_labels);
[locs,anatomy] = reconcile_locs_ieeg(clean_labels,clean_loc_labels,locs,anatomy);

%% locs
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

if ~isequal(labels,pc_labels)
    error('label mismatch');
end

%% Make inv dist matrix
A = inverse_dist(out_locs);
A = wrap_or_unwrap_adjacency_fc_toolbox(A);
A(nan_chs,:) = [];
A(:,nan_chs) = [];
labels(nan_chs) = [];


%%
figure
set(gcf,'position',[10 10 1000 500])
tiledlayout(1,2,'tilespacing','tight','padding','tight')

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

% re-vectorize
Avec = wrap_or_unwrap_adjacency_fc_toolbox(A);
pcvec = wrap_or_unwrap_adjacency_fc_toolbox(pc_net);
r = corr(Avec,pcvec,'rows','pairwise')

end