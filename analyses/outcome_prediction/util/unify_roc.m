function [unify_x,unify_y] = unify_roc(all_roc)

nx = 10000;
unify_x = linspace(0,1,nx);
unify_y = nan(length(all_roc),nx);

for i = 1:length(all_roc)
    curr_roc = all_roc{i};
    curr_x = curr_roc(:,1);
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