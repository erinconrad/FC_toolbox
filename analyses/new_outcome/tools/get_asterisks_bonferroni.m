function ast = get_asterisks_bonferroni(p,corr)

if p < 0.001/corr
    ast = '***';
elseif p < 0.01/corr
    ast = '**';
elseif p < 0.05/corr
    ast = '*';
else
    ast = '';
end

end