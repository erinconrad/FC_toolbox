function T = spike_atlas

%% Parameters
do_leader = 1;
which_atlas = 'brainnetome';% %'aal';'aal_bernabei';
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
if do_leader
    spikes = out.spike_leader_atlas;
else
    spikes = out.spikes_atlas;
end


%% Get soz loc-lat combos
soz_lats = out.all_soz_lats;
soz_locs = out.all_soz_locs;
soz_lat_loc = cellfun(@(x,y) [x,' ',y],soz_lats,soz_locs,'UniformOutput',false);

%% Remove patients without sz localizations
no_loc = cellfun(@isempty,soz_locs);
spikes(:,no_loc) = [];
soz_lats(no_loc) = [];
soz_locs(no_loc) = [];
soz_lat_loc(no_loc) = [];
npts = size(spikes,2);

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

%% Table with number of pts in each category
n_in_cat = n_things_in_cat(soz_broad,soz_regions);
T = array2table(n_in_cat','VariableNames',soz_regions);

%% Show spikes perc by patient
if 0
    turn_nans_gray(perc_spikes_broad)
    xticks(1:size(perc_spikes_broad,2))
    xticklabels(soz_broad)
    yticks(1:size(perc_spikes_broad,1))
    yticklabels(broad_regions)
end

%% SOZ loc -> spike loc
same_loc = soz_loc_predict_spike_loc(soz_broad,broad_locs,perc_spikes_broad,broad_regions);
same_lat = soz_lat_predict_spike_lat(soz_broad,broad_lats,perc_spikes_broad,broad_regions);

%% Pretty text
fprintf(['Spikes were most frequent in the hemisphere of the SOZ in %d out of %d (%1.1f%%) patients.'...
    ' This was similar for both left-sided seizure localization (%d of %d (%1.1f%%) patients) '...
    'and right-sided seizure localization (%d of %d (%1.1f%%) patients).\n'],...
    sum(same_lat(:,1)),sum(same_lat(:,2)),sum(same_lat(:,1))/sum(same_lat(:,2))*100,...
    same_lat(1,1),same_lat(1,2),same_lat(1,1)/same_lat(1,2)*100,...
    same_lat(2,1),same_lat(2,2),same_lat(2,1)/same_lat(2,2)*100);

%% Pretty text
fprintf(['Spikes were most frequent in the localization (mesial temporal vs temporal neocortical vs other cortex) of the SOZ in %d out of %d (%1.1f%%) patients.'...
    ' %d of %d (%1.1f%%) patients with mesial temporal lobe seizures had concordant spikes, '...
    '%d of %d (%1.1f%%) patients with temporal neocortical seizures had concordant spikes, and '...
    '%d of %d (%1.1f%%) patients with seizures in other cortex had concordant spikes.\n'],...
    sum(same_loc(:,1)),sum(same_loc(:,2)),sum(same_loc(:,1))/sum(same_loc(:,2))*100,...
    same_loc(1,1),same_loc(1,2),same_loc(1,1)/same_loc(1,2)*100,...
    same_loc(2,1),same_loc(2,2),same_loc(2,1)/same_loc(2,2)*100,...
    same_loc(3,1),same_loc(3,2),same_loc(3,1)/same_loc(3,2)*100);

%% Spikes loc -> SOZ loc
sz_same_loc = spike_loc_predict_soz_loc(broad_locs,broad_regions,perc_spikes_broad,soz_broad);
sz_same_lat =spike_lat_predict_soz_lat(broad_lats,broad_regions,perc_spikes_broad,soz_broad);

%% Pretty text

fprintf(['\nWe also asked how often the SOZ was concordant with spikes. ']);

fprintf(['%d out of %d (%1.1f%%) patients with left-predominant spikes had unilateral left sided seizure onset.'...
    ' (%d of %d (%1.1f%%) patients with right-predominant spikes had unilateral right sided seizure onset. \n'],...
    sz_same_lat(1,1),sz_same_lat(1,2),sz_same_lat(1,1)/sz_same_lat(1,2)*100,...
    sz_same_lat(2,1),sz_same_lat(2,2),sz_same_lat(2,1)/sz_same_lat(2,2)*100);

%% Pretty text
fprintf(['%d of %d (%1.1f%%) patients with mesial temporal lobe-predominant spikes had mesial temporal-only seizure onset, '...
    '%d of %d (%1.1f%% patients with temporal neocortical-predominant spikes had concordant seizure onset, and '...
    '%d of %d (%1.1f%%) patients with other cortex-predominant spikes had concordant seizure onset.\n'],...
   sz_same_loc(1,1),sz_same_loc(1,2),sz_same_loc(1,1)/sz_same_loc(1,2)*100,...
    sz_same_loc(2,1),sz_same_loc(2,2),sz_same_loc(2,1)/sz_same_loc(2,2)*100,...
    sz_same_loc(3,1),sz_same_loc(3,2),sz_same_loc(3,1)/sz_same_loc(3,2)*100);



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

%% ML
thing = perc_spikes_broad';
z = (thing - nanmean(thing,1))./nanstd(thing,[],1);
pt_idx = repmat((1:npts)',1,nb);
soz_bin = nan(npts,nb);
region_idx = repmat(1:6,npts,1);

for ip = 1:npts
    curr_soz = soz_broad{ip};
    soz_bin(ip,:) = strcmp(broad_regions,curr_soz);
end

unclear = sum(soz_bin,2) == 0;

% remove any unclear
thing(unclear,:) = [];
z(unclear,:) = [];
pt_idx(unclear,:) = [];
soz_bin(unclear,:) = [];
region_idx(unclear,:) = [];

% vectorize
thing_vec = thing(:);
pt_vec = pt_idx(:);
soz_vec = soz_bin(:);
z_vec = z(:);
region_vec = region_idx(:);

T = table(soz_vec,pt_vec,thing_vec,z_vec,region_vec);
T.pt_vec = nominal(T.pt_vec);
T.region_vec = (nominal(T.region_vec));

glme1 = fitglme(T,'soz_vec ~ thing_vec + (1|pt_vec)',...
    'Distribution','Poisson','Link','log');


glme2 = fitglme(T,'soz_vec ~ thing_vec + (1|pt_vec) + (1|region_vec)',...
    'Distribution','Poisson','Link','log');

%% Alt classifier
% Predict 1 of 9 SOZ localizations using spike pattern
thing = perc_spikes_broad';
class = soz_broad;
loc = soz_locs;
lat = soz_lats;
lat(~strcmp(lat,'left') & ~strcmp(lat,'right')) = {'bilateral'};
loc(contains(loc,'temporal')) = {'temporal'};
loc(~contains(loc,'temporal')) = {'other'};
%loc(~strcmp(loc,'mesial temporal') & ~strcmp(loc,'temporal neocortical')) = {'other'};


all_nan = sum(~isnan(thing),2) == 0;
thing(all_nan,:) = [];
class(all_nan) = [];
loc(all_nan) = [];
lat(all_nan) = [];
thing(isnan(thing)) = 0;

T = array2table(thing);
T = addvars(T,class);
T = addvars(T,lat);
T = addvars(T,loc);

return

%% Raw cluster without PCA
%{
perc_spikes_broad_pca = perc_spikes_broad;
perc_spikes_broad_pca(isnan(perc_spikes_broad_pca)) = 0;
thing = perc_spikes_broad_pca';

% get SSE as a function of k
SSE = find_optimal_k(thing,1:10);

% elbow plot
figure
plot(SSE,'-')
%}

%% PCA
perc_spikes_broad_pca = perc_spikes_broad;
perc_spikes_broad_pca(isnan(perc_spikes_broad_pca)) = 0;
[coeff,score,latent,tsquared,explained,mu] = pca(perc_spikes_broad_pca');

% first 3 components explain 80% of the variance in the data

%% Do K means clustering on first three components
thing = score(:,1:3);

% get SSE as a function of k
[SSE,idx_k] = find_optimal_k(thing,1:10);

% elbow plot
%plot(SSE,'-')

% natural elbows are 3 or 5. Try 5 since it's closer to what I want
k = 3;
idx = idx_k{k};
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
if 0
figure
turn_nans_gray(ksz)
xticks(1:nsozs)
xticklabels(soz_regions)
yticks(1:k)
end


    
end