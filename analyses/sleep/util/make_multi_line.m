function out = make_multi_line(text)

out = text;

for i = 1:length(text)
    
    curr_text = text{i};
    new_text = strrep(curr_text,' ',char(10));
    out{i} = new_text;
end


end