function lateralize_epilepsy

%% Parameters
do_plots = 1;
which_atlas = 'aal_bernabei';%'aal_bernabei';%%'brainnetome';
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

%% First, build symmetric coverage atlas
[symm_cov_atlas,all_bilateral] = build_symmetric_coverage_atlas(atlas,locs,lats);

% Also get symmetric spike maps
symm_spikes = nan(size(spikes));
for ip = 1:npts
    curr_spikes = spikes(:,ip);
    curr_bilateral = all_bilateral(:,ip);
    curr_spikes(~curr_bilateral) = nan;
    symm_spikes(:,ip) = curr_spikes;
end

%% Build a SOZ - non SOZ laterality ordered atlas
[soz_non_soz_ordered_atlas,soz,non_soz] = build_soz_ordered_atlas(atlas,left,right,right_lat,left_lat);
soz_non_soz_ordered_atlas_symm_cov = build_soz_ordered_atlas(symm_cov_atlas,left,right,right_lat,left_lat);

% same for spikes and bin_soz
soz_symm_spikes = nan(size(spikes));
soz_symm_bin_soz = nan(size(bin_soz));

for ip = 1:npts
    curr_spikes = symm_spikes(:,ip);
    curr_bin_soz = bin_soz(:,ip);
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    
    if curr_soz_right == 1
        soz_order = [find(right);find(left)];
    elseif curr_soz_left == 1
        soz_order = [find(left);find(right)];
    end
    
    soz_symm_spikes(1:length(soz_order),ip) = curr_spikes(soz_order);
    soz_symm_bin_soz(1:length(soz_order),ip) = curr_bin_soz(soz_order);
end

%% Get intrinsic connectivity of soz side and non soz side
lu_square = soz_non_soz_ordered_atlas_symm_cov(soz,soz,:);
rl_square = soz_non_soz_ordered_atlas_symm_cov(non_soz,non_soz,:);
lu_mean = squeeze(nanmean(lu_square,[1 2]));
rl_mean = squeeze(nanmean(rl_square,[1 2]));
soz_non_soz = [lu_mean rl_mean];    

% Average spike rate of soz and not
soz_spike_avg = nanmean(soz_symm_spikes(soz,:),1);
non_soz_spike_avg = nanmean(soz_symm_spikes(non_soz,:),1);
spikes_soz_non = [soz_spike_avg' non_soz_spike_avg'];

%% Show my SOZ-non SOZ ordered matrices
if do_plots
    figure
    set(gcf,'position',[10 10 1000 800])
    tiledlayout(3,2,'tilespacing','compact','padding','tight')

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
    plot_paired_data(soz_non_soz',{'SOZ side','non-SOZ side','non-SOZ side'},'Intra-hemispheric connectivity','paired',plot_type);
    title('Intrinisic connectivity in SOZ vs non-SOZ')
    assert(isequal(isnan(soz_non_soz(:,1)),isnan(soz_non_soz(:,2))))
    num_non_soz_higher = sum(diff(soz_non_soz,[],2)>0);
    num_non_soz_lower_or_equal = sum(diff(soz_non_soz,[],2)<=0);
    perc_non_soz_higher = 100*num_non_soz_higher/(num_non_soz_higher+num_non_soz_lower_or_equal);
    
    print(gcf,[plot_folder,'lat_',which_atlas],'-dpng')
    
end


%%

end