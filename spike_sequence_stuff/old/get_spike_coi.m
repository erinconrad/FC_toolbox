function [coi_ch,global_coi] = get_spike_coi(gdf,nchs,fs)

t2_seconds = 50e-3;
t2 = t2_seconds * fs;


coi_cell = cell(nchs,1);
for i = 1:nchs
    coi_cell{i} = [];
end
global_coi = nan(size(gdf,1),1);

for s = 1:size(gdf,1)
    time = gdf(s,2);
    ch = gdf(s,1);
    time_diff = abs(time - gdf(:,2));
    
    % find those less than t2 and not same ch
    close_enough = time_diff < t2 & gdf(:,1) ~= ch;
    
    % fill
    coi_cell{ch} = [coi_cell{ch};sum(close_enough)];
    global_coi(s) = sum(close_enough);
end

% average across all spikes for each ch
coi_ch = cellfun(@mean,coi_cell);
global_coi = mean(global_coi);

end