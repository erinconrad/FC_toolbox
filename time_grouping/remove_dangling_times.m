function [data,all_adj_bad] = remove_dangling_times(data,nruns)

 % find bad periods
all_bad = find(sum(~isnan(data),1) == 0);

% get periods adjacent to bad periods
before_bad = all_bad - 1;
after_bad = all_bad + 1;
all_adj_bad = unique([before_bad,after_bad]);

% remove those outside allowable times
all_adj_bad(all_adj_bad==0) = [];
all_adj_bad(all_adj_bad == nruns+1) = [];

% make these nans as well
data(:,all_adj_bad) = nan;


end