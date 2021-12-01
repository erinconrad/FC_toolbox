function mcout = test_ranking(rates,soz,nb)

npts = length(soz);
no_soz = cellfun(@isempty,soz);

%% Get true median ranking
all_rankings = nan(npts,1);
for ip = 1:npts
    
    % assign nan rates to be zero
    curr_rates = rates{ip};
    curr_rates(isnan(curr_rates)) = 0;
    
    % sort the rates in descending order
    [curr_rates,I] = sort(curr_rates,'descend');
    ranks = 1:length(curr_rates);
    ranks(I) = ranks;
    
    % get the ranks of the soz elecs
    curr_soz = soz{ip};
    soz_ranks = median(ranks(curr_soz));
    all_rankings(ip) = soz_ranks;
end

median_ranking_true = median(all_rankings(~no_soz));

%% Get MC median rankings
median_ranking_mc = nan(nb,1);
for ib = 1:nb
    if mod(ib,100) == 0
        fprintf('\nDoing %d out of %d\n',ib,nb);
    end
    temp_all_rankings = nan(npts,1);
    
    for ip = 1:npts
    
        % assign nan rates to be zero
        curr_rates = rates{ip};
        curr_rates(isnan(curr_rates)) = 0;

        % sort the rates in descending order
        [curr_rates,I] = sort(curr_rates,'descend');
        ranks = 1:length(curr_rates);
        ranks(I) = ranks;

        % get the number of soz electrodes
        curr_soz = soz{ip};
        nsoz = length(curr_soz);
        nelecs = length(curr_rates);
        
        % choose a random sample (without replacement) of electrodes equal
        % in number to nsoz
        fake_soz = randsample(nelecs,nsoz);
        
        soz_ranks = median(ranks(fake_soz));
        temp_all_rankings(ip) = soz_ranks;
    end
    median_ranking_mc(ib) = median(temp_all_rankings(~no_soz));
    
end

pval = (sum(median_ranking_mc <= median_ranking_true)+1)/(nb+1);

if 0
    plot(sort(median_ranking_mc),'o')
    hold on
    plot(xlim,[median_ranking_true median_ranking_true],'linewidth',2)
    title(get_p_text(pval))
    legend({'Monte Carlo','True'})
end

mcout.pval = pval;
mcout.all_rankings = all_rankings;
mcout.median_ranking_mc = median_ranking_mc;
mcout.median_ranking_true = median_ranking_true;


end