function [N,chi2] = generate_bootstrap_distribution(X)

totalN = sum(X,'all');
p_row = sum(X, 2)/sum(X(:));
p_col = sum(X, 1)/sum(X(:));
N = zeros(size(X,1),size(X,2));
for count = 1:totalN
   i = weighted_rand(size(X,1),p_row);
   j = weighted_rand(size(X,2),p_col);
   N(i,j) = N(i,j) + 1;


end

% calculate chi2
e = sum(N,2)*sum(N)/sum(N(:));
chi2 = (N-e).^2./e;
chi2 = nansum(chi2(:));

end

function b = weighted_rand(a,weights)

b = randsample(a,1,true,weights);

end