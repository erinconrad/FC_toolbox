function indices = reconcile_labels(all_labels)

n_labels = length(all_labels);

%% initialize cell of indices to keep
indices = cell(n_labels,1);

%% Get the labels common to all
common_labels = all_labels{1};
for i = 2:n_labels
    common_labels = intersect(common_labels,all_labels{i});
end

%% Go through and find the matching indices for each label set
for i = 1:n_labels
    
    % find the indices of the current set corresponding to the common set
    [~,locb] = ismember(common_labels,all_labels{i});
    
    % Make sure I can rederive common labels with this index
    if ~isequal(common_labels,all_labels{i}(locb)), error('what'); end
    
    indices{i} = locb;
end

end