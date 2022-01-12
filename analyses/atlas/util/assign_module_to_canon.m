function canon_assignments = assign_module_to_canon(avg_atlas,gamma,aal_to_yeo)

nreg = size(avg_atlas,1);
ncanon = size(aal_to_yeo,2);

% get communities
Ci=modularity_und(avg_atlas,gamma);
ncomms = length(unique(Ci));

% initialize list of canonical assignments for each region
canon_assignments = zeros(nreg,ncanon);

% loop over communities
for i = 1:ncomms
    
    % build vector with 1s for regions in that community
    x = zeros(nreg,1);
    x(Ci==i) = 1;
    
    % initialize overlap vector between x and each canonical system
    overlap = nan(ncomms,1);
    
    % Loop over canonical systems
    for j = 1:ncanon
        overlap(j) = corr(x,aal_to_yeo(:,j));
    end
    
    % find the canonical system with max overlap
    [~,I] = max(overlap);
    
    % assign the regions in that community to that canonical system
    canon_assignments(logical(x),I) = 1;
    
end

if 0
    imagesc(canon_assignments)
end


end