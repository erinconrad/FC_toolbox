function [all_ranks,all_soz_ranks,successes] = simple_rate_rank(rates,sozs)

%% Initialize data things
npts = length(rates);
all_ranks = cell(npts,1);
all_soz_ranks = cell(npts,1);
nchance = nan(npts,1);
all = nan(npts,1);
successes = nan(npts,1);

for i = 1:npts
    curr_soz = sozs{i};
    curr_rates = rates{i};
    
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
    soz_ranks = ranks(curr_soz);
    successes(i) = ranking_binomial_test(ranks,curr_soz);
    all_ranks{i} = ranks;
    all_soz_ranks{i} = soz_ranks;

end