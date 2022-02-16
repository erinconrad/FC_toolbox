function [x,pmf] = generalized_binomial_distribution(p)

pmf = 1;
n = length(p);

for i = 1:n
    pmf = conv(pmf,[1-p(i);p(i)]);
end

x = (0:n)';

end