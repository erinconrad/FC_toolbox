function [pval,nchance] = general_binomial_test(input)

% This function takes an nx2 "input" vector, where n is the number of
% independent tests, the first column is a 1 or a 0 indicating success or
% failure, and the second column is the probability of success for the
% test. It returns a p value for the question of how often a result as
% significant would be seen under the null hypothesis

%{


%}


%% Generate null distribution
% probabilities of getting each possible number of success under the null
% hypothesis
[x,pmf] = generalized_binomial_distribution(input(:,2));

%% True number of successes
%n = size(input,1);
n_failures = sum(input(:,1) == 0);


%% Two-tailed test (Wikipedia version)

this_pmf = pmf(x == n_failures);
pval = sum(pmf(pmf <= this_pmf));

%% Number successes expected by chance
nchance = sum(input(:,2));


end