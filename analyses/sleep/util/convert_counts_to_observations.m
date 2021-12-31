function observations = convert_counts_to_observations(counts,edges)

observations = nan(sum(counts),1);
ind = 0;

for i = 1:length(counts)
    N = counts(i);
    observations(ind+1:ind+N) = edges(i);
    ind = ind+N; 
end

end

%{

function observations = convert_counts_to_observations(counts,edges)

observations = nan(sum(counts),1);
ind = 0;

for i = 1:length(counts)
    N = counts(i);
    observations(ind+1:ind+N) = edges(i);
    ind = ind+N; 
end

end

%}