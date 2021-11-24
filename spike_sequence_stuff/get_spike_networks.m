function [dist,leader_ana] = get_spike_networks(summ,disc)

%% To do
%{ 
- think about optimizing gamma (right now using 1)
- decide the localization of the cluster (mode location??)
- see what they do with sleep
- get actual dot product

- find lower and higher RL channels in each cluster, see if their spike
rates go up with sleep. If lower goes up, then that means more independent
spikes. If higher goes up disproportionately, then more spread.
- correlation between RL in each cluster and spike increase with sleep

%}

%% Parameters
range = 0.2:0.01:1.8;
min_spikes = 1e3;

%% Add BCT toolbox
locations = fc_toolbox_locs;
addpath(genpath(locations.bct));
addpath(genpath(locations.script_folder))

%% Get basic info
coa = summ.coa;
locs = summ.locs;
rl = summ.rl;
locs = summ.locs;
ad = summ.ad;
ana = summ.ana_loc;
spikes = summ.spikes;
labels = summ.labels;

%% unwrap to 3 dimensional matrix
coa = wrap_or_unwrap_adjacency_fc_toolbox(coa);

%% Find and remove non-intracranial
ekg = find_non_intracranial(labels);
coa = coa(~ekg,~ekg,:);
locs = locs(~ekg,:);
rl = rl(~ekg,:);
ana = ana(~ekg);
spikes = spikes(~ekg,:);
labels = labels(~ekg);
ad = ad(~ekg,:);

ad = nanmean(ad,1);
avg_spikes = nanmean(spikes,2);

mean_rl = nanmean(rl,2);

 %% Determine "wake" and "sleep" times
% normalized ad
ad_norm = (ad - nanmedian(ad))./iqr(ad);
wake = ad_norm > disc;
sleep = ad_norm <= disc;


% sum over all times (nan because times that I throw out have no spikes)
sum_coa = nansum(coa,3);

% find optimal gamma
%find_optimal_gamma(sum_coa,range);

% cluster the network
[M,Q]=community_louvain(sum_coa);

% Get spikey communities and get distribution over time
nclusters = length(unique(M));
nchs = size(coa,1);
nruns = size(coa,3);
total_spikes = nan(nclusters,1);
dist = nan(nclusters,nruns);
leader_ana = cell(nclusters,1);
leader_chs = nan(nclusters,2);
rl_inc_corr = nan(nclusters,1);
for c = 1:nclusters
    sub_net = zeros(nchs,nchs);
    sub_net(M == c,M == c) = 1;
    sub_net = logical(sub_net);
    total_spikes(c) = sum(sum(sum_coa(sub_net)))/sum(M==c); % divide by number of elecs
    
    if 0
        imagesc(sub_net)
        pause
    end
    
    for r = 1:nruns
        % Get current network
        curr_net = coa(:,:,r);
        
        % take dot product of this network and current network
        dot_prod = sum(sum(curr_net .* sub_net));
        dist(c,r) = dot_prod;
    end
    
    % find the electrode with the lowest recruitment latency in that
    % cluster
    [~,lowest_rl] = min(mean_rl(M==c));
    
    % find the electrode with the most spikes in that cluster
    [~,spikiest] = max(avg_spikes(M==c));
    
    % get anatomical locations
    curr_ana = ana(M==c);
    leader_ana{c} = curr_ana{spikiest};
    
    chs = 1:nchs;
    curr_chs = chs(M==c);
    leader_chs(c,:) = [curr_chs(spikiest),curr_chs(lowest_rl)];
    
    % RL within cluster
    rl_in_cluster = mean_rl(M==c);
    
    % Spikes in cluster
    spikes_in_cluster = spikes(M==c,:);
    
    spikes_ws = [nanmean(spikes_in_cluster(:,wake),2),nanmean(spikes_in_cluster(:,sleep),2)];
    spikes_change = (spikes_ws(:,2)-spikes_ws(:,1))./spikes_ws(:,1);
    
    % correlation between RL and rate increase in sleep - if spikes occur
    % later in the cluster, do they have a bigger or smaller rate increase
    % with sleep?
    r = corr(rl_in_cluster,spikes_change,'rows','pairwise','type','spearman');
    rl_inc_corr(c) = r;
    
end

% Only care about spikey ones
spikey = total_spikes >= min_spikes;
nspikey = sum(spikey);

% function show them
if 0
    figure
    tiledlayout(nspikey,3)
    
    for c = 1:nclusters
        if spikey(c) == 0
            continue
        end
        
        sub_net = zeros(nchs,nchs);
        sub_net(M == c,M == c) = 1;
        
        nexttile
        imagesc(sub_net)
        
        nexttile
        plot(dist(c,:))
        
        nexttile
        
        elecs = M ==c;
        scatter3(locs(elecs,1),locs(elecs,2),locs(elecs,3),100,'ro','filled')
        hold on
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'ko')
        title(leader_ana{c})
        
    end
    
end
    
%% For output    
dist = dist(spikey,:);
leader_ana = leader_ana(spikey);

end