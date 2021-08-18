function stim_chs = return_stim_chs(cceps)

stim_chs = zeros(length(cceps.chLabels),1);
for i = 1:length(cceps.elecs)
    if ~isempty(cceps.elecs(i).arts)
        stim_chs(i) = 1;
    end
end

stim_chs = logical(stim_chs);

end