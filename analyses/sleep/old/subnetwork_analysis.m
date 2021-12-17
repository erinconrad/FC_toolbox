function out = subnetwork_analysis(disc)

%% Parameters
main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
nlocs = length(main_locs);
main{1} = main_locs;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

all_clust_ws = nan(nlocs,npts,2);
all_inter_clust_corr = nan(npts,1);

%% Loop over patients
for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    ad = summ.ad;
    labels = summ.labels;
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    
    %% Get spike networks
    [dist,leader_ana] = get_spike_networks(summ,disc);
    nclusters = size(dist,1);
    
     %% Determine "wake" and "sleep" times
    % normalized ad
    ad_norm = (ad - nanmedian(ad))./iqr(ad);
    wake = ad_norm > disc;
    sleep = ad_norm <= disc;
    
    %% Get average sleep vs wake distribution for each of the subnetworks
    % Also get the inter-cluster correlation
    if nclusters == 0, continue; end
    clust_ws = nan(nclusters,2);
    
    if nclusters > 1
        inter_clust_corr = nan(nchoosek(nclusters,2),1);
    end
    count = 0;
    for c = 1:nclusters
        dist_sleep = squeeze(nanmean(dist(c,sleep)));
        dist_wake = squeeze(nanmean(dist(c,wake)));
        clust_ws(c,:) = [dist_wake,dist_sleep];
        
        if nclusters > 1
            for d = 1:c-1
                count = count + 1;
                inter_clust_corr(count) = corr(dist(c,:)',dist(d,:)','rows','pairwise');
            end
        end
    end
    
    % get average inter-cluster correlation for the patient and fill array
    if nclusters > 1
        all_inter_clust_corr(p) = nanmean(inter_clust_corr);
    end
    
    %% Fill up distributions based on locs
    for l = 1:nlocs
        % find the clusters matching this loc
        ic = strcmp(leader_ana,main_locs{l});
        
        % average distributions for these clusters
        curr_avg = nanmean(clust_ws(ic,:),1);
        
        % fill up
        all_clust_ws(l,p,:) = curr_avg;
        
    end

end

out.all_clust_ws = all_clust_ws;
out.main = main;

if 0
    plot_type = 'errorbar';
   interaction_plot_and_stats(all_clust_ws,make_multi_line(main{1}),'Distribution',...
    {'Awake','Asleep'},0,plot_type);
title('Wake/sleep spike rate by anatomical location') 
end

if 0
figure
nclust = size(all_clust_ws,1);
for l = 1:nlocs
    for i = 1:2
        plot(l+randn(nclust,1)*0.01+0.1,squeeze(all_clust_ws(:,l,1)),'o')
        hold on
        plot(l+randn(nclust,1)*0.01-0.1,squeeze(all_clust_ws(:,l,2)),'o')
    end
end
end


end