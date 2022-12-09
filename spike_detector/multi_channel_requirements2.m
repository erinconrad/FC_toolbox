function final_spikes = multi_channel_requirements2(gdf,nchs,fs)

%% Parameters
min_chs = 2; % spike should be on at least 2 channels
max_chs = nchs * 0.5; % on no more than half the channels
min_time = 100 * 1e-3 * fs; % 100 ms to look for other spikes

final_spikes = [];

s = 1;
curr_seq = s;
last_time = gdf(s,2);

while s<size(gdf,1)
    
    % move to next spike time
    new_time = gdf(s+1,2);
    
    % if it's within the time diff
    if new_time - last_time < min_time
        curr_seq = [curr_seq;s+1]; % append it to the current sequence
        
        if s == size(gdf,1) - 1
            
            % done with sequence, check if the number of involved chs is
            % appropriate
            %l = length(curr_seq);
            l = length(unique(gdf(curr_seq,1)));
            if l >= min_chs && l <= max_chs
                final_spikes = [final_spikes;gdf(curr_seq,:) ...
                    gdf(curr_seq,2)-min(gdf(curr_seq,2))];
            end
            
        end
    else
        % done with sequence, check if the length of sequence is
        % appropriate
        %l = length(curr_seq);
        l = length(unique(gdf(curr_seq,1)));
        if l >= min_chs && l <= max_chs
            final_spikes = [final_spikes;gdf(curr_seq,:) gdf(curr_seq,2)-min(gdf(curr_seq,2))];
        end
        
        % reset sequence
        curr_seq = s+1;
    end
    
    % increase the last time
    last_time = gdf(s+1,2);
    
    % increase the current spike
    s = s+1;
    
end



end