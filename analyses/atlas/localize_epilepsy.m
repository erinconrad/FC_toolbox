function localize_epilepsy

%% Parameters
which_atlas = 'aal_bernabei'; %'aal';
plot_type = 'scatter';
broad_regions = {'left mesial temporal','right mesial temporal',...
    'left temporal neocortical','right temporal neocortical',...
    'left other cortex','right other cortex'};
nb = length(broad_regions);

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

if strcmp(which_atlas,'aal_bernabei'), names = names'; end

%% Get soz atlas identifiers
sozs = out.sozs;

%% Get soz broader localizations
soz_out = load('out.mat');
soz_out = soz_out.out.circ_out;
soz_locs = soz_out.all_locs;
soz_lats = soz_out.all_lats;
soz_lat_loc = cellfun(@(x,y) [x,' ',y],soz_lats,soz_locs,'UniformOutput',false);
if 0
    table(soz_locs,soz_lats)
end

if 0
%% Get electrodes to ignore
ignore = regions_to_ignore(which_atlas,names);

atlas(ignore,:,:) = nan;
atlas(:,ignore,:) = nan;
end

%% Localize regions into broad categories
broad = localize_regions(names,which_atlas);
if 0
    table(names,broad)
end

%% Ignore other
%{
ignore = cellfun(@isempty,broad);
atlas(ignore,:,:) = nan;
atlas(:,ignore,:) = nan;
%}

%% Get the average intrinsic connectivity of each of the broad regions and get soz
broad_connectivity = nan(nb,npts);
soz_broad = zeros(nb,npts);
for ib = 1:nb
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
zbroad = (broad_connectivity - nanmean(broad_connectivity,2))./nanstd(broad_connectivity,[],2);

if 0
    figure
    nexttile
    turn_nans_gray(zbroad)
    
    nexttile
    turn_nans_gray(soz_broad)
end

%% For each patient, compare the regional connectivity of SOZ region to non-SOZ regions
soz_not = nan(npts,2);
soz_not_non_normalized = nan(npts,2);
for ip = 1:npts
    curr_soz_broad = soz_broad(:,ip);
    
    if sum(curr_soz_broad) == 0
        continue;
    end
    soz_not(ip,1) = (zbroad(curr_soz_broad,ip));
    soz_not(ip,2) = nanmean(zbroad(~curr_soz_broad,ip));
    
    soz_not_non_normalized(ip,:) = [broad_connectivity(curr_soz_broad,ip),...
        nanmean(broad_connectivity(~curr_soz_broad,ip))];
end


figure
set(gcf,'position',[221 425 1220 372])
tiledlayout(1,2)
nexttile
stats = plot_paired_data(soz_not_non_normalized',{'SOZ','non-SOZ','non-SOZ'},'Raw intrinsic connectivity','paired',plot_type);

nexttile
stats = plot_paired_data(soz_not',{'SOZ','non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);




%% Test if on an atlas-region level, the NS is different for SOZ than others
% Get node strengths (mean across edges for each node).
% I do mean instead of sum because nansum of all nans is 0, but nanmean of
% all nans is nan
if 0
ns = squeeze(nanmean(atlas,1));

% Convert ns to z-score, comparing each region to that of the same region in other patients
z = (ns - nanmean(ns,2))./nanstd(ns,[],2);

% For each patient, get the average z score of the SOZ electrodes and non-SOZ electrodes

z_soz_no = nan(npts,2);
for ip = 1:npts
    curr_z = z(:,ip);
    curr_sozs = sozs{ip};
    is_soz = ismember(atlas_nums,curr_sozs);
    z_soz = nanmean(curr_z(is_soz));
    z_not = nanmean(curr_z(~is_soz));
    z_soz_no(ip,:) = [z_soz z_not];
    
    
    if 0
        fprintf('\n%s SOZ:\n',pt_names{ip});
        names(is_soz)
        pause
        
    end
    
end

% Remove rows with any nans
any_nans = any(isnan(z_soz_no),2);
z_soz_no(any_nans,:) = [];

stats = plot_paired_data(z_soz_no',{'SOZ','Non-SOZ','non-SOZ'},'Node strength (z-score)','paired',plot_type);

fprintf(['\nThe mean edge weights of the SOZ regions (median %1.2f) was '...
    'lower than that of the non-SOZ regions (median %1.2f) (T+ = %1.2f, p = %1.3f).\n'],...
    stats.medians(1),stats.medians(2),stats.Tpos,stats.pval);

fprintf(['\nThe connectivity was higher in the non-SOZ for %d (%1.1f%%) of patients.\n'],...
    sum(diff(z_soz_no,1,2)>0),sum(diff(z_soz_no,1,2)>0)/size(z_soz_no,1)*100);
end


end