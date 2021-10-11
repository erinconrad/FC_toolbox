function [coa,rl] = build_sequences(gdf,nchs,fs)

%% How close spikes need to be
t2_seconds = 50e-3;
t2 = t2_seconds * fs;

%% Initialize array
coa = zeros(nchs,nchs);
rl = cell(nchs,1);

%% Loop over spikes
for s = 1:size(gdf,1)
    time = gdf(s,2);
    ch = gdf(s,1);
    
    % difference in spike time between this one and all others
    time_diff = abs(time - gdf(:,2));
    
    % find those less than t2 and not same ch. I will say these spikes are
    % in the same spike sequence
    close_enough = time_diff < t2 & gdf(:,1) ~= ch;
 
    % Get the sequence start time and the latency of this channel in the
    % sequence
    seq_start_time = min([gdf(close_enough,2);time]);
    seq_ch_time = (time - seq_start_time)/fs;
    
    % find the channels that are in the sequence with it
    close_chs = gdf(close_enough,1);
    close_chs = unique(close_chs);
    close_chs(close_chs == ch) = [];
    
    % Add count to coa
    coa(ch,close_chs) = coa(ch,close_chs) + 1;
    
    % add time to rl if at least two channels
    if sum(close_enough) > 0
        rl{ch} = [rl{ch};seq_ch_time];
    end

end

%% confirm coa symmetric
%assert(issymmetric(coa))
coa = min(coa,coa');
assert(issymmetric(coa));


%% Average across spikes for rl
rl = cellfun(@mean,rl);