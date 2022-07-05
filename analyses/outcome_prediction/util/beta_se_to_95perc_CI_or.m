function [or,or_ci] = beta_se_to_95perc_CI_or(beta,se)

or = exp(beta);
ci = [beta-1.96*se,beta+1.96*se];
or_ci = [exp(ci(1)),exp(ci(2))];

end