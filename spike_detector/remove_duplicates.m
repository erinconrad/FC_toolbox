function [gdf,n_removed] = remove_duplicates(gdf)
% max time diff
off = 50e-3;

% sort by time
gdf = sortrows(gdf,2);

keep = ones(size(gdf,1),1);

% take diff of times
diff_times = [inf;diff(gdf(:,2))];

% take diff of chs
diff_chs = [inf;diff(gdf(:,1))];

% find those that are close in time and the same ch
too_close = abs(diff_times) < off & diff_chs == 0;

keep(too_close) = 0;
keep = logical(keep);

n_removed = sum(~keep);
gdf(~keep,:) = [];


end