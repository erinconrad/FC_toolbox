function coa = build_coa_from_seq(seq)

nchs = size(seq,2);
nspikes = size(seq,1);
coa = zeros(nchs,nchs);
for s = 1:nspikes
    non_zero = find(seq(s,:) ~= 0);
    for i = 1:length(non_zero)
        for j = 1:i-1
            
            coa(non_zero(i),non_zero(j)) = coa(non_zero(i),non_zero(j)) + 1;
            coa(non_zero(j),non_zero(i)) = coa(non_zero(j),non_zero(i)) + 1;
        end
    end
end


end