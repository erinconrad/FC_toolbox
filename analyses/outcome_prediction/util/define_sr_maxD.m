function sr = define_sr_maxD(locs)

% Get interelectrode distance matrix
D = make_interdist_matrix(locs);

% Take maximum of D
maxD = max(D,[],'all');

% Define this to be the search radius
sr = maxD;

end
