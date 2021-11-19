function counts = sz_counts_bins(szs,nbins)

counts = zeros(nbins,1);
for z = 1:size(szs,1)
    counts(szs(z,2)) = counts(szs(z,2)) + 1;
end

end