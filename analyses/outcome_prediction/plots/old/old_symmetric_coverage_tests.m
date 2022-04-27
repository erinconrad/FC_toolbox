function old_symmetric_coverage_tests

%% Parameters
do_plots = 0;
which_atlas = 'brainnetome';%;'aal_bernabei';%%'brainnetome';
plot_type = 'scatter';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];

bct_folder= locations.bct;
out_folder = [results_folder,'analysis/atlas/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end
if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));


%% Load atlas and get region names and spikes
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;
atlas = out.atlas;
names = out.atlas_names;
spikes = out.spikes_atlas;
nregions = length(names);
assert(nregions==size(atlas,1))
npts = size(atlas,3);
sozs = out.sozs;
atlas_nums = out.atlas_nums;
bin_soz = (cell2mat(cellfun(@(x) ismember(atlas_nums',x),sozs,'uniformoutput',false)))';

%% Load soz lats
soz_lats = out.all_soz_lats;
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');

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
[symm_cov_atlas,all_bilateral] = build_symmetric_coverage_atlas(atlas,locs,lats);
atlas = symm_cov_atlas;

%% Double check symmetric coverage
% Looks good!
for ip = 1:npts
    curr = atlas(:,:,ip);
    avg_rows = nanmean(curr,2);
    
    % find non nan
    non_nan = find(~isnan(avg_rows));
    
    % confirm contralateral is non nan
    for i = 1:length(non_nan)
        assert(~isnan(avg_rows(contra_index(non_nan(i)))))
    end
    
end

%% Get average connectivity of SOZ (to other things) and that of contralateral region
soz_intra = nan(npts,2);
lr_soz_intra = nan(npts,2);
soz_all = nan(npts,2);
lr_soz_all = nan(npts,2);
hemi = nan(npts,2);
hemi_lr = nan(npts,2);
for ip = 1:npts
    
    %% SOZ connectivity
    % get soz indices
    curr_soz = bin_soz(:,ip);
    
    % get the regions contalateral to soz
    contra_soz = contra_index(curr_soz);
    contra_soz(isnan(contra_soz)) = [];
    bin_contra_soz = zeros(length(curr_soz),1);
    bin_contra_soz(contra_soz) = 1;
    bin_contra_soz = logical(bin_contra_soz);
    
    % confirm that difference in number of soz and contra soz is just those
    % regions that have no contralateral thing 
    switch which_atlas
        case 'aal_bernabei'
            assert(abs(sum(curr_soz)-sum(bin_contra_soz)) == ...
                sum(ismember([names(curr_soz)';names(bin_contra_soz)'],names(end-8:end))))
        case 'brainnetome'
            assert((sum(curr_soz)==sum(bin_contra_soz)))
    end
    
    
    %% Get the left and right
    all_ipsi_and_contra_soz = [find(curr_soz);find(bin_contra_soz)];
    left_soz = all_ipsi_and_contra_soz(left(all_ipsi_and_contra_soz) == 1);
    right_soz = all_ipsi_and_contra_soz(right(all_ipsi_and_contra_soz) == 1);
    
    assert(sum(strcmp(lats(left_soz),'L')) == length(left_soz) && ...
        sum(strcmp(lats(right_soz),'R'))==length(right_soz));
    
    lr_soz_intra(ip,:) =  [nanmean(atlas(left_soz,left_soz,ip),'all'),...
        nanmean(atlas(right_soz,right_soz,ip),'all')];
    lr_soz_all(ip,:) =  [nanmean(atlas(left_soz,:,ip),'all'),...
        nanmean(atlas(right_soz,:,ip),'all')];
    
    %% get SOZ-SOZ connectivity and contra-SOZ - contra-SOZ connectivity
    soz_intra(ip,:) = [nanmean(atlas(curr_soz,curr_soz,ip),'all'),...
        nanmean(atlas(bin_contra_soz,bin_contra_soz,ip),'all')];
        
    soz_all(ip,:) = [nanmean(atlas(curr_soz,:,ip),'all'),...
        nanmean(atlas(bin_contra_soz,:,ip),'all')];
    
    %% Holohemispjheric
    % get everything on the side of the soz
    if right_lat(ip) == 1
        ipsi_lats = strcmp(lats,'R');
        contra_lats = strcmp(lats,'L');
    elseif left_lat(ip) == 1
        ipsi_lats = strcmp(lats,'L');
        contra_lats = strcmp(lats,'R');
    else
        continue;
    end
    
    % Get soz side - soz side connectivity and contralateral connectivity
    hemi(ip,:) = [nanmean(atlas(ipsi_lats,ipsi_lats,ip),'all'),...
        nanmean(atlas(contra_lats,contra_lats,ip),'all')];
   
    hemi_lr(ip,:) = [nanmean(atlas(left,left,ip),'all'),...
        nanmean(atlas(right,right,ip),'all')];
end

%% Build a SOZ - non SOZ laterality ordered atlas
soz_non_soz_ordered_atlas = build_soz_ordered_atlas(atlas,left,right,right_lat,left_lat);

%% get confusion matrix
hemi_diff = hemi_lr(:,1) - hemi_lr(:,2);
predicted = cell(length(soz_lats),1);
predicted(hemi_diff > 0) = {'right'};
predicted(hemi_diff < 0) = {'left'};
empty = hemi_diff == 0 | isnan(hemi_diff) | (~strcmp(soz_lats,'right') & ~strcmp(soz_lats,'left'));
predicted(empty) = [];
lats_for_conf = soz_lats;
lats_for_conf(empty) = [];

if 0
    table(lats_for_conf,predicted)
end

conf_out = confusion_matrix(predicted,lats_for_conf,0);

%% COnfusion matrix to lateralize epilepsy
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

print(gcf,[plot_folder,'symm_pred_',which_atlas],'-dpng')

if do_plots

figure
set(gcf,'position',[-2023         710        1781         891])
tiledlayout(2,6,'tilespacing','compact','padding','tight')

nexttile([1 3])
pretty_matrix(all_bilateral(~neither_lat,:),...
    {'SOZ','non-SOZ'},sum(left),[],1);
colormap(gca,[0.5,0.5,0.5;1 1 1])
title('Regions with symmetric coverage')
xlabel('Patient')
ylabel('Region')

nexttile([1 3])
pretty_matrix(nanmean(soz_non_soz_ordered_atlas(~neither_lat,~neither_lat,:),3),...
    {'SOZ','non-SOZ'},sum(left),'r^2',0)
title('Average connectivity (symmetric coverage only)')
xlabel('Patient')

nexttile([1 2])
plot_paired_data(hemi',{'SOZ side','non-SOZ side','non-SOZ side'},'Intra-hemispheric connectivity','paired',plot_type);
title('Intra-hemispheric connectivity on side of SOZ vs non-SOZ')

nexttile([1 2])
plot_paired_data(soz_all',{'SOZ','contralateral region','contralateral region'},'Average connectivity','paired',plot_type);
title('Average connectivity to SOZ vs contralateral region')

nexttile([1 2])
plot_paired_data(soz_intra',{'SOZ','contralateral region','contralateral region'},'Intrinsic connectivity','paired',plot_type);
title('Intrinsic connectivity in SOZ vs contralateral region')

print(gcf,[plot_folder,'symm_',which_atlas],'-dpng')


%% LR to check
figure
set(gcf,'position',[-2023         710        1781         350])
tiledlayout(1,3,'tilespacing','compact','padding','tight')
nexttile
plot_paired_data(hemi_lr',{'left','right','right'},'Intra-hemispheric connectivity','paired',plot_type);
title('Intra-hemispheric connectivity on left vs right')

nexttile
plot_paired_data(lr_soz_all',{'left','right','right'},'Average connectivity','paired',plot_type);
title('Average connectivity to SOZ localization (left vs right)')

nexttile
plot_paired_data(lr_soz_intra',{'left','right','right'},'Intrinsic connectivity','paired',plot_type);
title('Intrinsic connectivity in SOZ localization (left vs right)')

print(gcf,[plot_folder,'symm_lr_',which_atlas],'-dpng')



end


end