function fine_localization

%{
At individual region level, get node strength. Then normalize across
patients. Do SOZ regions have higher or lower connectivity than non SOZ
regions?

Edit this to:
- normalize at level of EDGES
- see if average intrinsic connectivity of SOZ laterality different from
that of non-SOZ
- then restrict to correct laterality.
- see if node strength different
%} 

%% Parameters
which_atlas = 'brainnetome';%'aal_bernabei';%'brainnetome';% %'aal';'aal_bernabei';
plot_type = 'scatter';
coverage_limit = 10;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];
bct_folder= locations.bct;
out_folder = [results_folder,'analysis/atlas/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load atlas
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;

%% Load soz lats
soz_out = load('out.mat');
soz_out = soz_out.out.circ_out;
soz_lats = soz_out.all_lats;
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');

atlas = out.atlas;
names = out.atlas_names;
pt_names = out.pt_names;
atlas_nums = out.atlas_nums;
nregions = length(names);
assert(nregions==size(atlas,1))
npts = size(atlas,3);
nelecs = out.n_elecs_all;
sozs = out.sozs;
bin_soz = (cell2mat(cellfun(@(x) ismember(atlas_nums',x),sozs,'uniformoutput',false)))';


%% Indicate regions with low coverage
any_coverage = squeeze(any(~isnan(atlas),2));

% find regions with low coverage
low_coverage = sum(any_coverage,2) < coverage_limit;


if 0
    figure
    nexttile
    turn_nans_gray(any_coverage)
    
    nexttile
    turn_nans_gray(any_coverage(~low_coverage,:))
    yticks(1:sum(~low_coverage))
    yticklabels(names(~low_coverage))
end

%% Restrict atlas and soz
%{
atlas = atlas(~low_coverage,~low_coverage,:);
bin_soz = bin_soz(~low_coverage,:);
names = names(~low_coverage);

%}
%% Remove white matter
if 1
white_matter = ismember(names,'White_Matter');
names = names(~white_matter);
bin_soz = bin_soz(~white_matter,:);
atlas = atlas(~white_matter,~white_matter,:);
end

if 0
    figure
    nexttile
    turn_nans_gray(nanmean(atlas,3))
    
    nexttile
    turn_nans_gray(bin_soz)
end

%% Get locs and lats for atlas names
[locs,lats,loc_nums] = lateralize_regions(names,which_atlas);

%% Normalize atlas
atlas = (atlas - nanmean(atlas,3))./nanstd(atlas,[],3);

%% Come up with new order so that I can visualize network according to SOZ/non-SOZ
% First, order by left and then right
new_order = reorder_lr(locs,lats);
reordered_lats = lats(new_order);
reordered_locs = locs(new_order);
new_order_left = strcmp(reordered_lats,'L');
new_order_right = strcmp(reordered_lats,'R');
new_order_neither = ~new_order_left & ~new_order_right;
reordered_atlas = atlas(new_order,new_order,:);
lat_labels = cell(246,1);
lat_labels(1:123) = {'SOZ'};
lat_labels(124:end) = {'Non'};
loc_labels = reordered_locs;

soz_order_atlas = nan(size(atlas));
for ip = 1:npts
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    
    if curr_soz_right == 0 && curr_soz_left == 0, continue; end
    
    if curr_soz_left == 1
        soz_order_atlas(:,:,ip) = reordered_atlas(:,:,ip);
        
    elseif curr_soz_right == 1
        soz_order_atlas(1:sum(new_order_right),1:sum(new_order_right),ip) = reordered_atlas(new_order_right,new_order_right,ip);
        soz_order_atlas(sum(new_order_right)+1:sum(new_order_right)+sum(new_order_left),...
            sum(new_order_right)+1:sum(new_order_right)+sum(new_order_left),ip) = reordered_atlas(new_order_left,new_order_left,ip);
        soz_order_atlas(sum(new_order_right)+sum(new_order_left)+1:end,...
            sum(new_order_right)+sum(new_order_left)+1:end,ip) = reordered_atlas(new_order_neither,new_order_neither,ip);
        
    end
    
end


if 1
    turn_nans_gray(nanmean(soz_order_atlas,3))
    yticks(1:size(soz_order_atlas,1))
    yticklabels(cellfun(@(x,y) [x,' ',y],lat_labels,loc_labels,'uniformoutput',false))
end

%% Ns
ns = squeeze(nanmean(atlas,2));
z = ns;
%z = (ns-nanmean(ns,2))./nanstd(ns,[],2);

if 0
    figure
    nexttile
    turn_nans_gray(z)
    yticks(1:size(z,1))
    yticklabels(names)
    
    nexttile
    turn_nans_gray(bin_soz)
end

%% Plot orders
%plot_orders_mats(ns,bin_soz);

%% Z soz vs not
soz_not = nan(npts,2);
soz_not_raw = nan(npts,2);
for ip = 1:npts
    curr_soz = bin_soz(:,ip);
    soz_not(ip,:) = [nanmean((z(curr_soz,ip))), nanmean((z(~curr_soz,ip)))];
    soz_not_raw(ip,:) = [nanmean((ns(curr_soz,ip))), nanmean((ns(~curr_soz,ip)))];
end

if 1
figure
nexttile
stats = plot_paired_data((soz_not)',{'SOZ','non-SOZ','non-SOZ'},'Normalized connectivity','paired',plot_type);

nexttile

stats = plot_paired_data((soz_not_raw)',{'SOZ','non-SOZ','non-SOZ'},'Raw connectivity','paired',plot_type);
end



end