function [pval,exp] = binomial_generalized(ps,true)
    pmf = [1];
    for i = 1:length(ps)
        pmf = conv(pmf, [1-ps(i); ps(i)]);
    end
    x = (0:length(ps))';
    
    index = x >= true;
    pval = sum(pmf(index)); % shoudl this be 2 sided instead???
    exp = mean(ps);
end