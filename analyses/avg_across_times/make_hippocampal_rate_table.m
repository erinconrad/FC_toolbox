function make_hippocampal_rate_table(out)

all_left = [];
all_right = [];
names = {};

for i = 1:length(out)
    left_hipp_rate = nanmean(out(i).avg_spikes(out(i).left_right_hipp(:,1)==1));
    right_hipp_rate = nanmean(out(i).avg_spikes(out(i).left_right_hipp(:,2)==1));
    
    all_left = [all_left;left_hipp_rate];
    all_right = [all_right;right_hipp_rate];
    names = [names;out(i).name];
    
end

avg_left = nanmean(all_left);
avg_right = nanmean(all_right);
max_left = max(all_left);
max_right = max(all_right);

T = table(names,all_left,all_right);

writetable(T,'jim_table.csv');

end