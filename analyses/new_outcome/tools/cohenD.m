function d = cohenD(pop1,pop2)
% remove nans
pop1(isnan(pop1)) = [];
pop2(isnan(pop2)) = [];

n1 = length(pop1);
n2 = length(pop2);
v1 = var(pop1);
v2 = var(pop2);

s = sqrt(((n1-1)*v1+(n2-1)*v2)/(n1+n2-2));
d = (mean(pop1)-mean(pop2))/s;

end