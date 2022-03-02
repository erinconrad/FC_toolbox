function spike_atlas

%% Parameters
which_atlas = 'brainnetome';%'aal_bernabei';%'brainnetome';% %'aal';'aal_bernabei';
plot_type = 'scatter';

broad_locs = {'mesial temporal','temporal neocortical','other cortex'};
broad_lats = {'left','right'};


%% Make all broad regions from locs and lats
nlocs = length(broad_locs);
nlats = length(broad_lats);
broad_regions = cell(nlocs*nlats,1);
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

%% Get atlas info
names = out.atlas_names;
spikes = out.spikes_atlas;
npts = size(spikes,2);

%% Get soz loc-lat combos
soz_lats = out.all_soz_lats;
soz_locs = out.all_soz_locs;
soz_lat_loc = cellfun(@(x,y) [x,' ',y],soz_lats,soz_locs,'UniformOutput',false);

% Reassign into one of 9 categories
[soz_broad,soz_regions] = assign_soz_localization(soz_lat_loc,broad_regions);
nsozs = length(soz_regions);

%% Localize regions into broad categories
broad = localize_regions(names,which_atlas);

%% Average spikes in each category
spikes_broad = nan(nb,npts);
for ib = 1:nb
    
    % get atlas regions in that broad region
    rs = strcmp(broad,broad_regions{ib});
    
    % average spike rates in those atlas regions
    avg_spikes = nanmean(spikes(rs,:),1);
    
    % assign to new matrix
    spikes_broad(ib,:) = avg_spikes;
    
end

%% Also get percent of patient's spikes in that region
perc_spikes_broad = spikes_broad./nansum(spikes_broad,1)*100;

%% group patients across 9 groups and remake matrices
soz_spikes = nan(nb,nsozs);
soz_spikes_perc = nan(nb,nsozs);
for isoz = 1:nsozs
    curr_pts = strcmp(soz_broad,soz_regions{isoz});
    soz_spikes(:,isoz) = nanmean(spikes_broad(:,curr_pts),2);
    soz_spikes_perc(:,isoz) = nanmean(perc_spikes_broad(:,curr_pts),2);
end

%% Show
if 0
figure
set(gcf,'position',[-1076 741 1800 600])
tiledlayout(1,2)

nexttile
turn_nans_gray(soz_spikes)
xticks(1:nsozs)
xticklabels(soz_regions)
yticks(1:nb)
yticklabels(broad_regions)
title('Average spike rate')
set(gca,'fontsize',15)
c = colorbar;
ylabel(c,'spikes/elecs/min','fontsize',15)
xlabel('SOZ localization')
ylabel('Spike region')


nexttile
turn_nans_gray(soz_spikes_perc)
xticks(1:nsozs)
xticklabels(soz_regions)
yticks(1:nb)
yticklabels(broad_regions)
title('% Spikes in region')
set(gca,'fontsize',15)
c = colorbar;
ylabel(c,'%','fontsize',15)
xlabel('SOZ localization')
ylabel('Spike region')
end

%% PCA
perc_spikes_broad_pca = perc_spikes_broad;
perc_spikes_broad_pca(isnan(perc_spikes_broad_pca)) = 0;
[coeff,score,latent,tsquared,explained,mu] = pca(perc_spikes_broad_pca');

% first 3 components explain 80% of the variance in the data

%% Do K means clustering on first three components
thing = score(:,1:3);

% get SSE as a function of k
SSE = find_optimal_k(thing);

% elbow plot
%plot(SSE,'-')

% natural elbows are 3 or 5. Try 5 since it's closer to what I want
k = 5;
idx = kmeans(thing,k);
if 0
figure
for i = 1:k
    plot3(thing(idx==i,1),thing(idx==i,2),thing(idx==i,3),'o','linewidth',2);
    hold on
end
end

%% Look at cluster assignments for each category
ksz = nan(k,nsozs);
for ic = 1:k
    for is = 1:nsozs
        ksz(ic,is) = sum(idx==ic & strcmp(soz_broad,soz_regions{is}));
    end
end

%% Plot cluster assignment by SOZ category
figure
turn_nans_gray(ksz)
xticks(1:nsozs)
xticklabels(soz_regions)
yticks(1:k)

    
end