function n_in_cat = n_things_in_cat(things,cats)

n_in_cat = nan(length(cats),1);

for i = 1:length(cats)
    curr = cats{i};
    n_in_cat(i) = sum(strcmp(things,curr));
    
end

end