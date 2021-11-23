function [dist,M,spikey] = get_spike_networks(summ)

%% To do
%{ 
- think about optimizing gamma (right now using 1)
- decide the localization of the cluster (mode location??)
- see what they do with sleep
- get actual dot product

%}

%% Parameters
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
ana = summ.ana_loc;

mean_rl = nanmean(rl,2);

% unwrap to 3 dimensional matrix
coa = wrap_or_unwrap_adjacency_fc_toolbox(coa);

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
for c = 1:nclusters
    sub_net = zeros(nchs,nchs);
    sub_net(M == c,M == c) = 1;
    sub_net = logical(sub_net);
    total_spikes(c) = sum(sum(sum_coa(sub_net)));
    
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
    
    % get anatomical locations
    curr_ana = ana(M==c);
    leader_ana{c} = curr_ana{lowest_rl};
    
end

% Only care about spikey ones
spikey = total_spikes >= min_spikes;
nspikey = sum(spikey);

% function show them
if 1
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
    
    


end