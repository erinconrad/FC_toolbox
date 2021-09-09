function overlap = find_overlapping_groups(pairs)

pair1 = pairs(1,:);
pair2 = pairs(2,:);

% put them in order!!!!!!

% if pair1 entirely before pair2
if pair1(2) <= pair2(1)
    overlap = 0;
% if pair1 entirely before pair2
elseif pair2(2) <= pair1(1)
    overlap = 0;
else
    overlap = 1;
end


end