function regional_localization_epilepsy

%% Parameters
which_atlas = 'brainnetome';%'aal_bernabei';%'brainnetome';% %'aal';'aal_bernabei';
plot_type = 'scatter';

broad_locs = {'mesial temporal','temporal neocortical','other cortex'};
broad_lats = {'left','right'};
broad_regions = cell(6,1);
nb = length(broad_regions);
count = 0;
for i = 1:length(broad_locs)
    for j = 1:length(broad_lats)
        count = count + 1;
        broad_regions{count} = [broad_lats{j},' ',broad_locs{i}];
    end
end
same_side = nan(nb,2);
for ib = 1:nb
    C  = strsplit(broad_regions{ib},' ');
    lat = C{1};
    curr_same_side = find(contains(broad_regions,lat));
    curr_same_side(curr_same_side == ib) = [];
    same_side(ib,:) = curr_same_side;
end



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

atlas = out.atlas;
names = out.atlas_names;
pt_names = out.pt_names;
atlas_nums = out.atlas_nums;
nregions = length(names);
assert(nregions==size(atlas,1))
npts = size(atlas,3);

%% Get soz loc-lat combos
soz_lats = out.all_soz_lats;
soz_locs = out.all_soz_locs;
soz_lat_loc = cellfun(@(x,y) [x,' ',y],soz_lats,soz_locs,'UniformOutput',false);


%% Localize regions into broad categories
broad = localize_regions(names,which_atlas);

%% Get the average intrinsic connectivity of each of the broad regions and get soz
broad_connectivity = nan(nb,npts);
soz_broad = zeros(nb,npts);
% Loop over broad regions (6 total)
for ib = 1:nb
    % Loop over patients
    for ip = 1:npts
        
        % Get the regions in that broad region
        curr_broad = broad_regions{ib};
        in_region = strcmp(broad,curr_broad);
        
        % get intrinsic connectivity of that region
        intrinsic = nanmean(atlas(in_region,in_region,ip),'all');
        broad_connectivity(ib,ip) = intrinsic;
        
        % Get soz
        curr_soz = soz_lat_loc{ip};
        soz_broad(ib,ip) = strcmp(curr_soz,curr_broad);
        
    end
end
assert(sum(sum(soz_broad,1)<=1) == length(soz_broad))
soz_broad = logical(soz_broad);

%% Get normalized regional connectivity
% For each patient and each region, I ask how intrinsically connected is
% that region relative to the same region in other patients.
zbroad = (broad_connectivity - nanmean(broad_connectivity,2))./nanstd(broad_connectivity,[],2);

%% For each patient, compare the normalized regional connectivity of SOZ region to non-SOZ regions ON SAME SIDE
soz_not = nan(npts,2);
soz_not_non_normalized = nan(npts,2);
for ip = 1:npts
    curr_soz_broad = soz_broad(:,ip); % get this patient's SOZ region
    
    if sum(curr_soz_broad) == 0 % skip if no SOZ (half the patients!)
        continue;
    end
    
    curr_same_side = same_side(curr_soz_broad,:);
    
    soz_not(ip,1) = (zbroad(curr_soz_broad,ip)); % get normalized connectivity of SOZ
    soz_not(ip,2) = nanmean(zbroad(curr_same_side,ip)); % normalized connectivity of other regions, averaged across regions
    
    soz_not_non_normalized(ip,:) = [broad_connectivity(curr_soz_broad,ip),...
        nanmean(broad_connectivity(curr_same_side,ip))]; %same but not normalized connectivity (not controlling for normal anatomical differences)
    
   
end

figure
set(gcf,'position',[221 425 1220 372])
tiledlayout(1,2)
nexttile
stats = plot_paired_data(soz_not_non_normalized',{'SOZ','non-SOZ','non-SOZ'},'Raw intrinsic connectivity','paired',plot_type);
title('Raw regional connectivity by SOZ status')

nexttile
stats = plot_paired_data(soz_not',{'SOZ','non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);
title('Normalized regional connectivity by SOZ status')

end