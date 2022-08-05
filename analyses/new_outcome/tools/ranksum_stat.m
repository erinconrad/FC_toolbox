function [p,W] = ranksum_stat(continuous_thing,binary_thing)

if sum(binary_thing) == 0
    W = nan;
    p = nan;
    return
end

if sum(~isnan(continuous_thing)) == 0
    W = nan;
    p = nan;
    return
end
[p,~,stats] = ranksum(continuous_thing(binary_thing==1),continuous_thing(binary_thing==0));
W = stats.ranksum;

end