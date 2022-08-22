function sub_bin = subsample_chs(chs_bin,n)

chs = find(chs_bin);
sub_bin = zeros(length(chs_bin),1);

% Return if n == 0
if n == 0
    return;
end

if length(chs) > 1 % if more than 1, subsample appropriately
    sub = randsample(chs,n);
    sub_bin(sub) = 1;
else
    
    sub_bin = chs_bin; % if only 1, keep it
    assert(sum(chs_bin)==1)
    
end
   
end