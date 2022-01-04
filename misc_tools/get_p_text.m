function ptext = get_p_text(p)

if p < 0.001
    ptext = '\itp\rm < 0.001';
elseif p < 0.05
    ptext = sprintf('\\itp\\rm = %1.3f',p);
else
    ptext = sprintf('\\itp\\rm = %1.2f',p);
end

end