function conc = concordance_at_the_top(vec1,vec2)

assert(length(vec1)==length(vec2))
n = length(vec1);
conc = nan(n,1);

for i = 1:n
    conc(i) = length(intersect(vec1(1:i),vec2(1:i)))/i;
end

end