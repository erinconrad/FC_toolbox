function ptext = get_p_html(p)

if p < 0.001
    ptext = '<i>p</i> < 0.001';
elseif p < 0.05
    ptext = sprintf('<i>p</i> = %1.3f',p);
else
    ptext = sprintf('<i>p</i> = %1.2f',p);
end

end