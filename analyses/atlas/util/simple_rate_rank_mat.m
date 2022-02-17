function [nranks,median_soz_ranks,successes] = simple_rate_rank_mat(rates,sozs,successes)

%% Initialize data things
npts = size(rates,2);
median_soz_ranks = nan(npts,1);
nranks = nan(npts,1);
orig_pts = 1:npts;

soz_chance = nan(npts,2);

for i = 1:npts
    curr_soz = sozs(:,i);
    curr_rates = rates(:,i);
    
    % skip if nan
    if sum(curr_soz) == 0, continue; end
    
    % remove nans
    nan_things = isnan(curr_rates);
    curr_soz(nan_things) = [];
    curr_soz = logical(curr_soz);    
    curr_rates(nan_things) = [];
    
    % Get rank of SOZ
    % sort the rates in descending order
    [~,I] = sort(curr_rates,'descend');
    ranks = 1:length(curr_rates);
    ranks(I) = ranks;
    nranks(i) = length(ranks);
    soz_ranks = ranks(curr_soz);
    
    
    median_soz_ranks(i) = median(soz_ranks);
    
    soz_chance(i,:) = [median(soz_ranks) median(ranks)];

end


% remove nans
nan_rows = any(isnan(soz_chance),2);
orig_pts(nan_rows) = [];
nranks(nan_rows) = [];
median_soz_ranks(nan_rows) = [];
soz_chance(nan_rows,:) = [];
alt_successes = (soz_chance(:,1)>soz_chance(:,2));
alt_pval = 2*binocdf(sum(alt_successes==0),length(alt_successes),0.5);

%% Alt two sided p-value
% Test both directions
%{
nhigher = sum(soz_chance(:,1) > (soz_chance(:,2)));
npts = sum(~any(isnan(soz_chance),2));

% test higher
if nhigher > npts/2
    pval_binom = 2*binocdf(sum(soz_chance(:,1) <= sum(soz_chance(:,2))),npts,0.5);
% test lower
else
    pval_binom = 2*binocdf(sum(soz_chance(:,1) >= sum(soz_chance(:,2))),npts,0.5);
end
%}

end