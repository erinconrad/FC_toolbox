function sz_blocks = get_sz_blocks(times,sz_times)

sz_blocks = nan(size(sz_times,1),2);
for is = 1:size(sz_times,1)
    [~,I] = min(abs(times-sz_times(is,1)));
    if times(I) > sz_times(is,1), I = I-1; end
    sz_blocks(is,1) = I;
    
    [~,I] = min(abs(times-sz_times(is,2)));
    if times(I) < sz_times(is,2), I = I+1; end
    sz_blocks(is,2) = I;
    
end

end