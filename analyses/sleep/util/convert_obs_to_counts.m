function counts = convert_obs_to_counts(obs,nbins)

counts = zeros(nbins,1);

for i = 1:length(obs)
    counts(obs(i)) = counts(obs(i)) + 1;
end

end