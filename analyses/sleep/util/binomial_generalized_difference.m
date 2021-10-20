function pval = binomial_generalized_difference(ps,true_diff)
    pmf = [1];
    for i = 1:length(ps)
        pmf = conv(pmf, [1-ps(i); ps(i)]);
    end
    pmf = conv(pmf, flipud(pmf));
    x = (-length(ps):length(ps))';
    
    index = abs(x) >= abs(true_diff);
    pval = sum(pmf(index));
end