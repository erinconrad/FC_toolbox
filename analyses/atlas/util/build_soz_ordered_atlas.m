function [soz_ordered_atlas,soz,non_soz] = build_soz_ordered_atlas(lr_ordered_atlas,left,right,right_soz,left_soz)

npts = size(lr_ordered_atlas,3);
soz_ordered_atlas = nan(size(lr_ordered_atlas));

soz = zeros(size(lr_ordered_atlas,1),1);
non_soz = zeros(size(lr_ordered_atlas,1),1);
soz(left) = 1;
non_soz(right) = 1;

soz = logical(soz);
non_soz = logical(non_soz);

for ip = 1:npts
    curr_soz_right = right_soz(ip);
    curr_soz_left = left_soz(ip);
    curr_atlas = lr_ordered_atlas(:,:,ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    
    if curr_soz_right == 1
        soz_order = [find(right);find(left)];
    elseif curr_soz_left == 1
        soz_order = [find(left);find(right)];
    end
    
    soz_ordered_atlas(1:length(soz_order),1:length(soz_order),ip) = curr_atlas(soz_order,soz_order);
    
    
end

end