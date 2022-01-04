function out = ranking_binomial_test(ranks,soz)

median_soz_ranking = median(ranks(soz));
median_overall_ranking = median(ranks);

% it's a success if the median soz ranking is less than the median overall
% ranking (are soz electrodes spikier than chance?)
if median_soz_ranking < median_overall_ranking
    out = 1;
else
    out = 0;
end

end