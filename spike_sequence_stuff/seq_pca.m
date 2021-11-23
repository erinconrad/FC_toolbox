function seq_pca(seq)

% binarize: 1 if co-spike, 0 if not.
seq(~isnan(seq)) = 1;
seq(isnan(seq)) = 0;

% test: do standard coa
coa = build_coa_from_seq(seq);

% pca
[W,H,D] = nnmf(seq,5);




end