function out = ranking_binomial_test(ranks,soz)

median_soz_ranking = median(ranks(soz));
median_overall_ranking = median(ranks);

if median_soz_ranking < median_overall_ranking
    out = 1;
else
    out = 0;
end

end