function find_optimal_gamma(sum_coa,range)

Q = nan(length(range),1);

for ir = 1:length(range)
    gamma = range(ir);
    [~,Q(ir)]=community_louvain(sum_coa,gamma);
end

plot(range,Q)

end