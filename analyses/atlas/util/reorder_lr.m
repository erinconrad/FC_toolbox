function new_order = reorder_lr(locs,lats)

% reorders atlas regions from original order to left on top then right then
% whatever is left over 

% Do all left and then all right

nr = length(locs);
new_order = nan(nr,1);

% find left and right
left = strcmp(lats,'L');
right = strcmp(lats,'R');
neither = ~left & ~right;

% fill up left and then right and then remaining
new_order(1:sum(left)) = find(left);
new_order(sum(left)+1:sum(left)+sum(right)) = find(right);
new_order(sum(left)+sum(right)+1:end) = find(neither);

end