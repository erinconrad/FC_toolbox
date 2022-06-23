function [unify_x,unify_y] = unify_roc(all_roc)

%% Take a bunch of individual model ROCs and fit to a single ROC curve.
% Make a range of X values (10,000) total, find the closest X in the
% current model and fill in the y with the corresponding y

nx = 10000;
unify_x = linspace(0,1,nx);
unify_y = nan(length(all_roc),nx);

% Loop over models
for i = 1:length(all_roc)
    curr_roc = all_roc{i};
    if isempty(curr_roc), continue; end
    curr_x = curr_roc(:,1);
    
    % Loop over output (unified xs)
    for ix = 1:nx
        curr_ix = unify_x(ix);
        
        % find the curr x closest to this ix
        [~,I] = min(abs(curr_ix-curr_x));
        
        % fill the y in with the corresponding y
        corr_y = curr_roc(I,2);
        unify_y(i,ix) = corr_y;
        
    end
end

end