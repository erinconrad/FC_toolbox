function fine_localization

%{
At individual region level, get node strength. Then normalize across
patients. Do SOZ regions have higher or lower connectivity than non SOZ
regions?

Edit this to:
- see if average intrinsic connectivity of Left different from right - no,
same
- see if average intrinsic connectivity of SOZ laterality different from
that of non-SOZ (don't normalize)
- then restrict to correct laterality.
- then exclude regions without many patients to normalize
- then normalize
- see if node strength different for SOZ vs non SOZ locations
%} 

%% Parameters
which_atlas = 'aal_bernabei';'brainnetome';%'aal_bernabei';%'brainnetome';% %'aal';'aal_bernabei';
plot_type = 'scatter';
coverage_limit = 30;

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


%% Get locs and lats for atlas names
[locs,lats,loc_nums] = lateralize_regions(names,which_atlas);
left = strcmp(lats,'L');
right = strcmp(lats,'R');

%% Get left and right intrinsic connectivity
left_intrinsic = squeeze(nanmean(atlas(left,left,:),[1 2]));
right_intrinsic = squeeze(nanmean(atlas(right,right,:),[1 2]));

%% Sanity check (I expect no): Is left intrinsic connectivity diff from right?
% No
if 0
figure
stats = plot_paired_data(([left_intrinsic right_intrinsic])',{'left','right','right'},'Intrinsic connectivity','paired',plot_type);
title('Left vs right intrinsic connectivity')
end

%% Get SOZ and non-SOZ laterality intrinsic connectivity

soz_non_intrinsic = nan(npts,2);
for ip = 1:npts
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    
    if curr_soz_right == 1
        soz_non_intrinsic(ip,:) = [right_intrinsic(ip) left_intrinsic(ip)];
    elseif curr_soz_left == 1
        soz_non_intrinsic(ip,:) = [left_intrinsic(ip) right_intrinsic(ip)];
    end
end

%% I expect yes: Is SOZ intrinsic connectivity diff from non-SOZ?
% Yes
if 0
figure
stats = plot_paired_data(soz_non_intrinsic',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
title('SOZ vs non-SOZ intrinsic connectivity')
end

%% Confirmation of above test 1: for regions with symmetric coverage, is SOZ side intrinsic connectivity less?
all_bilateral = find_symmetric_coverage(atlas,lats,locs);
symm_soz_not = nan(npts,2);
for ip = 1:npts
    curr_bilateral = logical(all_bilateral(:,ip));
    if sum(curr_bilateral) == 0, continue; end
    curr_atlas = atlas(:,:,ip);
    curr_atlas(~curr_bilateral,:) = nan;
    curr_atlas(:,~curr_bilateral) = nan;
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    if curr_soz_right == 1
        symm_soz_not(ip,:) = [nanmean(curr_atlas(right,right),'all') nanmean(curr_atlas(left,left),'all')];
    elseif curr_soz_left == 1
        symm_soz_not(ip,:) = [nanmean(curr_atlas(left,left),'all') nanmean(curr_atlas(right,right),'all')];
    end
end

%% Confirmation of above test 2: is normalized SOZ intrinsic connectivity diff from non-SOZ?

if 1
figure
stats = plot_paired_data(symm_soz_not',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
title({'SOZ vs non-SOZ intrinsic connectivity','(Symmetric coverage only)'})
end

%% Remove regions with low coverage
% find regions with low coverage
low_coverage = sum(any_coverage,2) < coverage_limit;
atlas = atlas(~low_coverage,~low_coverage,:);
bin_soz = bin_soz(~low_coverage,:);
names = names(~low_coverage);

%% Normalize edges
z = (z-nanmean(z,3))./nanstd(z,[],3);

%% 

end