function circ_stats(data)

ngroups = size(data,1);
npts = size(data,2);

if npts == 2
    [~,p,~,stats] = ttest(data');
    
elseif npts > 2
    tbl = ranova(
else
    error('what');
end

end