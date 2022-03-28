function predict_y = predict_conn_from_loc(f,locs)

% make inter-distance matrix
D = make_interdist_matrix(locs);

% convert to 1D
D = wrap_or_unwrap_adjacency_fc_toolbox(D);
x = D;

predict_y = (f.p1 * x + f.p2)./(x + f.q1);

% wrap it back to 2D
predict_y = wrap_or_unwrap_adjacency_fc_toolbox(predict_y);

end