function [coa,rl,global_coi,seq_lengths,leader] = build_sequences(gdf,nchs,fs)

%% How close spikes need to be
t2_seconds = 50e-3;
t2 = t2_seconds * fs;

%% Initialize array
nspikes = size(gdf,1);
coa = zeros(nchs,nchs);
rl = cell(nchs,1);
global_coi = nan(size(gdf,1),1);
%seq_matrix = nan(nspikes,nchs);
seqs = {};

%% Loop over spikes
for s = 1:size(gdf,1)
    time = gdf(s,2);
    ch = gdf(s,1);
    
    % Take last sequence
    if ~isempty(seqs)
        last_seq = seqs{end};
        last_ch = last_seq(end,1);
        last_time = last_seq(end,2);
        if time < last_time
            error('what');
        end
        
        if time - last_time < t2
            seqs{end}(end+1,:) = [ch time];
        else
            seqs{end+1}(1,:) = [ch time];
        end
    else
        seqs{1}(1,:) = [ch time];
    end
    
    % difference in spike time between this one and all others
    time_diff = abs(time - gdf(:,2));
    %unsigned_time_diff = (time - gdf(:,2));
    
    % find those less than t2 and not same ch. I will say these spikes are
    % in the same spike sequence
    close_enough = time_diff < t2 & gdf(:,1) ~= ch;
    
    % How many other spikes co-occur with this spike?
    global_coi(s) = sum(close_enough);
 
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
    
    % Get the latency of all of the co-spiking channels for that spikes
    %{
    seq_matrix(s,gdf(close_enough,1)) = unsigned_time_diff(close_enough);
    seq_matrix(s,ch) = 0; % add itself
    %}
    
    % add time to rl if at least two channels
    if sum(close_enough) > 0
        rl{ch} = [rl{ch};seq_ch_time];
    end

end

seq_lengths = cellfun(@(x) size(x,1), seqs);

%% Get the leader electrode for each sequence
seq_leader = cellfun(@(x) x(1,1), seqs);


% Convert this to an nchx1 vector, mostly zeros, with number of times the
% electrode is the leader
leader = zeros(nchs,1);
for i = 1:length(seq_leader)
    leader(seq_leader(i)) = leader(seq_leader(i)) + 1;
end



%% confirm coa symmetric
%assert(issymmetric(coa))
coa = min(coa,coa');
assert(issymmetric(coa));


%% Average across spikes for rl
rl = cellfun(@mean,rl);
global_coi = mean(global_coi);