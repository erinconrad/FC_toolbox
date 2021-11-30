function out_text = add_asterisks_to_chosen(in_text,chosen)
    
out_text = in_text;
for i = 1:length(in_text)
    if chosen(i) == 1
        out_text{i} = [out_text{i},'***'];
    end
end


end