function pts_per_edge = measure_sparsity(atlas)

nedges = size(atlas,1);

pts_per_edge = nan(nedges,1);

for ie = 1:nedges
    edge = atlas(ie,:);
    pts_per_edge(ie) = sum(~isnan(edge));
end

end