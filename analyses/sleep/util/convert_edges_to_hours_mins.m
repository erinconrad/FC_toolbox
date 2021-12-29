function hours_mins = convert_edges_to_hours_mins(edges)

hours_mins = cell(length(edges),1);
for i = 1:length(edges)
    s = seconds(edges(i));
    [h,m] = hms(s);
    if h < 10
        htext = sprintf('0%d',h);
    else
        htext = sprintf('%d',h);
    end
    
    if m == 0
        mtext = '00';
    else
        mtext = sprintf('%d',m);
    end
    
    
    hours_mins{i} = sprintf('%s:%s',htext,mtext);
end