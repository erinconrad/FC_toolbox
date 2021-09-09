function ptext = get_p_text(p)

if p < 0.001
    ptext = 'p < 0.001';
elseif p < 0.05
    ptext = sprintf('p = %1.3f',p);
else
    ptext = sprintf('p = %1.2f',p);
end

end