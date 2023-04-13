function [out_labels,allowed_idx] = find_mt_symmetric_coverage(labels,allowable_labels)

%% Find possibly allowed labels (those that match allowable labels (mesial temporal-targeted electrodes)
allowed = ismember(labels,allowable_labels);
allowed_labels = labels(allowed);
nallowed = sum(allowed);

%% Of this restricted set, only allow those that have symmetric coverage
final_allowed_idx = zeros(nallowed,1);

for in = 1:nallowed
    % get current label
    curr = allowed_labels{in};

    % define contralateral one
    if contains(curr,'L')
        opp = strrep(curr,'L','R');
    elseif contains(curr,'R')
        opp = strrep(curr,'R','L');
    else
        error('waht');
    end

    % if contralateral one also there
    [ia,ib] = ismember(opp,allowed_labels);
    if ia == 1
        final_allowed_idx(in) = 1;
        final_allowed_idx(ib) = 1;
    end
end

final_allowed_idx = logical(final_allowed_idx);
out_labels = allowed_labels(final_allowed_idx);

allowed_idx = ismember(labels,out_labels);

if 0
    table(out_labels(contains(out_labels,'L')),out_labels(contains(out_labels,'R')))
end
% Make sure left and right ones are the same
assert(isequal(strrep(sort(out_labels(contains(out_labels,'L'))),'L',''),strrep(sort(out_labels(contains(out_labels,'R'))),'R','')))
  

end