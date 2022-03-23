function symmetric_coverage_tests

%% Parameters
do_plots = 1;
which_atlas = 'aal_bernabei';%%'brainnetome';
plot_type = 'scatter';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];

bct_folder= locations.bct;
out_folder = [results_folder,'analysis/atlas/'];
plot_folder = [results_folder,'analysis/atlas/plots/'];
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
soz_all = nan(npts,2);
hemi = nan(npts,2);
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
    
    % get SOZ-SOZ connectivity and contra-SOZ - contra-SOZ connectivity
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
   
    
end

%% Build a SOZ - non SOZ laterality ordered atlas
soz_non_soz_ordered_atlas = build_soz_ordered_atlas(atlas,left,right,right_lat,left_lat);

figure
set(gcf,'position',[-2023         710        1781         891])
tiledlayout(2,6,'tilespacing','compact','padding','tight')

nexttile([1 3])
pretty_matrix(all_bilateral(~neither_lat,:),...
    {'SOZ','non-SOZ'},sum(left),[],1)
title('Regions with symmetric coverage')

nexttile([1 3])
pretty_matrix(nanmean(soz_non_soz_ordered_atlas(~neither_lat,~neither_lat,:),3),...
    {'SOZ','non-SOZ'},sum(left),'r^2',0)
title('Average connectivity (symmetric coverage only)')

nexttile([1 2])
plot_paired_data(hemi',{'SOZ side','non-SOZ side','non-SOZ side'},'Intra-hemispheric connectivity','paired',plot_type);
title('Intra-hemispheric connectivity on side of SOZ vs non-SOZ')

nexttile([1 2])
plot_paired_data(soz_all',{'SOZ','contralateral region','contralateral region'},'Average connectivity','paired',plot_type);
title('Average connectivity to SOZ vs contralateral region')

nexttile([1 2])
plot_paired_data(soz_intra',{'SOZ','contralateral region','contralateral region'},'Intrinsic connectivity','paired',plot_type);
title('Intrinsic connectivity in SOZ vs contralateral region')

end