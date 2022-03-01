function lateralize_epilepsy

%% Parameters
do_plots = 1;
which_atlas = 'brainnetome';%'aal_bernabei';%%'brainnetome';
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


%% Load atlas and get region names
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;
atlas = out.atlas;
names = out.atlas_names;
nregions = length(names);
assert(nregions==size(atlas,1))

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

%% First, build symmetric coverage atlas
[symm_cov_atlas,all_bilateral] = build_symmetric_coverage_atlas(atlas,locs,lats);

%% Build a SOZ - non SOZ laterality ordered atlas
[soz_non_soz_ordered_atlas,soz,non_soz] = build_soz_ordered_atlas(atlas,left,right,right_lat,left_lat);
soz_non_soz_ordered_atlas_symm_cov = build_soz_ordered_atlas(symm_cov_atlas,left,right,right_lat,left_lat);


%% Get intrinsic connectivity of soz and non soz
lu_square = soz_non_soz_ordered_atlas_symm_cov(soz,soz,:);
rl_square = soz_non_soz_ordered_atlas_symm_cov(non_soz,non_soz,:);
lu_mean = squeeze(nanmean(lu_square,[1 2]));
rl_mean = squeeze(nanmean(rl_square,[1 2]));
soz_non_soz = [lu_mean rl_mean];    

%% Show my SOZ-non SOZ ordered matrices
if do_plots
    figure
    set(gcf,'position',[10 10 1000 800])
    tiledlayout(2,2,'tilespacing','compact','padding','tight')

    nexttile
    pretty_matrix(nanmean(soz_non_soz_ordered_atlas(~neither_lat,~neither_lat,:),3),...
        {'SOZ','non-SOZ'},sum(left),'r^2',0)
    title('Average connectivity')
    
    nexttile
    pretty_matrix(all_bilateral(~neither_lat,:),...
        {'SOZ','non-SOZ'},sum(left),[],1)
    title('Regions with symmetric coverage')
    
    
    nexttile
    pretty_matrix(nanmean(soz_non_soz_ordered_atlas_symm_cov(~neither_lat,~neither_lat,:),3),...
        {'SOZ','non-SOZ'},sum(left),'r^2',0)
    title('Average connectivity (symmetric coverage only)')


    nexttile
    stats = plot_paired_data(soz_non_soz',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
    title('Intrinisic connectivity in SOZ vs non-SOZ')
    
end

print(gcf,[plot_folder,'lat_',which_atlas],'-dpng')



end