function oops_symmetric_spike_tests

%% Parameters
do_plots = 1;
which_atlas = 'aal_bernabei';%%'brainnetome';
plot_type = 'scatter';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];

bct_folder= locations.bct;
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end
if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load out file and get roc stuff
out = load([out_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
rl = out.all_rl;
ns = out.all_ns;
soz = out.all_soz_bin;
soz_loc = out.all_soz_locs;
npts = length(soz);
labels = out.all_labels;
locs = out.all_locs;
fc = out.all_fc;
soz_lats = out.all_soz_lats;
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');


%% Load atlas and get region names and spikes
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

%% get atlas stuff
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;
spikes_atlas = atlas_out.spikes_atlas;
atlas_nums = atlas_out.atlas_nums;
atlas_names = atlas_out.atlas_names;
ns_atlas = atlas_out.atlas;
ns_atlas = squeeze(nanmean(ns_atlas,2));


%% Get locs and lats for atlas names
[locs,lats] = lateralize_regions(names,which_atlas);
left = strcmp(lats,'L');
right = strcmp(lats,'R');
neither_lat = ~left & ~right;

% confirm atlas has as many on right as left
assert(sum(left)==sum(right));

%% Re-order atlas to be left then right then neither
lr_order = reorder_lr(locs,lats);

left = left(lr_order);
right = right(lr_order);
neither_lat = neither_lat(lr_order);
atlas = atlas(lr_order,lr_order,:);
names = names(lr_order);
locs = locs(lr_order);
lats = lats(lr_order);
spikes = spikes(lr_order,:);
bin_soz = bin_soz(lr_order,:);
atlas_nums = atlas_nums(lr_order);

%% index of contralateral
contra_index = nan(size(left));
contra_index(1:sum(left)) = ([1:sum(left)])'+sum(left);
contra_index(sum(left)+1:sum(left)*2) = ([1:sum(left)])';

% make sure first half locs are same as last half locs
assert(isequal(locs(1:sum(left)),locs(sum(left)+1:sum(left)*2)))
assert(isequal(locs(contra_index(1:sum(left))),locs(contra_index(sum(left)+1:sum(left)*2))))

%% First, build symmetric coverage atlas
[~,all_bilateral] = build_symmetric_coverage_atlas(atlas,locs,lats);

%% First, build symmetric coverage atlas
% Also get symmetric spike maps
symm_spikes = nan(size(spikes));
for ip = 1:npts
    curr_spikes = spikes(:,ip);
    curr_bilateral = all_bilateral(:,ip);
    curr_spikes(~curr_bilateral) = nan;
    symm_spikes(:,ip) = curr_spikes;
end

%% SOZ vs non-SOZ hemispheric spike rate
lr_spikes = [nanmean(symm_spikes(left,:),1)',nanmean(symm_spikes(right,:),1)'];

%% get confusion matrix
hemi_diff = lr_spikes(:,1) - lr_spikes(:,2);
predicted = cell(length(soz_lats),1);
predicted(hemi_diff < 0) = {'right'}; % more spikes on the right -> right sided epilepsy
predicted(hemi_diff > 0) = {'left'};
empty = hemi_diff == 0 | isnan(hemi_diff) | (~strcmp(soz_lats,'right') & ~strcmp(soz_lats,'left'));
predicted(empty) = [];
lats_for_conf = soz_lats;
lats_for_conf(empty) = [];

if 0
    table(lats_for_conf,predicted)
end

conf_out = confusion_matrix(predicted,lats_for_conf,0);

figure
turn_nans_gray([1 0;0 1])
colormap(gca,[0.8500, 0.3250, 0.0980;0, 0.4470, 0.7410])
xticks(1:conf_out.nclasses)
xticklabels(conf_out.classes)
yticks(1:conf_out.nclasses)
yticklabels(conf_out.classes)
xlabel(conf_out.xlabel)
ylabel(conf_out.ylabel)
hold on
for i = 1:conf_out.nclasses
    for j = 1:conf_out.nclasses
        text(i,j,sprintf('%d',conf_out.mat(j,i)),'horizontalalignment','center','fontsize',25,'fontweight','bold')
    end
end
title(sprintf('Accuracy: %1.1f%%, PPV: %1.1f%%, NPV: %1.1f%%',conf_out.accuracy*100,...
    conf_out.ppv*100,conf_out.npv*100))
set(gca,'fontsize',15)

print(gcf,[plot_folder,'symm_pred_spikes_',which_atlas],'-dpng')

end