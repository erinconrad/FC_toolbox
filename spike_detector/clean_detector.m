function gdf = clean_detector(values,fs, options)

%{
Inputs:
- values: an nsamples x nchs array of eeg values
- fs: sampling rate in Hz

Outputs:
- gdf: an nspike x 5 array with data for each detected spike. Spikes are in
ascending order of time and columns are the following
   - column 1: spike channel
   - column 2: spike time in SAMPLES (not s)
   - column 3: spike duration in SAMPLES
   - column 4: spike amplitude in uV
   - column 5: time (in samples) of the spike relative to its spike train (first
   spike is 0)

Dependencies:
- eegfilt
- FindPeaks
- multi_channel_requirements2

Information:
- This was originally written by Camilo Bermudez 7/31/13 and modified by
Erin Conrad at UPenn 12/9/22.
- Updated arguments 3/9/2023, bscheid@seas.upenn.edu
%}

%% Parameters
arguments
    values
    fs
    options.tmul = 19; % minimum relative amplitude (compared to baseline)
    options.absthresh = 100; % minimum absolute amplitude (uV)
    options.sur_time = 0.5; % surround time (in s) against which to compare for relative amplitude
    options.close_to_edge = 0.05; % time (in s) surrounding start and end of sample to ignore 
    options.too_high_abs = 1e3; % amplitude above which I reject it as artifact
    options.spkdur = [15 200]; % spike duration must be within this range (in ms)
    options.lpf1 = 30; % low pass filter for spikey component
    options.hpf  = 7; % high pass filter for spikey component
end

    tmul = options.tmul; 
    absthresh = options.absthresh; 
    sur_time = options.sur_time; 
    close_to_edge = options.close_to_edge; 
    too_high_abs = options.too_high_abs; 
    spkdur = options.spkdur*fs/1000; 
    lpf1 = options.lpf1; 
    hpf  = options.hpf; 

%% Initialize things
all_spikes  = [];
nchs = size(values,2);

%% Iterate channels and detect spikes
for j = 1:nchs
    
    % initialize out array with final spike info
    out = [];
    
    % specify ch
    data = values(:,j);
    
    % Skip if all nans
    if sum(isnan(data)) > 0, continue; end
    
    % re-adjust the mean of the data to be zero
    data = data - nanmean(data);
        
    % initialize array with tentative spike info
    spikes   = [];

    % Low pass filter to remove artifact
    lpdata = eegfilt(data, lpf1, 'lp',fs); % low pass filter

    % high pass filter to get the spikey part
    hpdata   = eegfilt(lpdata, hpf, 'hp',fs); % high pass filter

    % establish the baseline for the relative amplitude threshold
    lthresh = median(abs(hpdata)); 
    thresh  = lthresh*tmul;     % this is the final threshold we want to impose

    % Run the spike detector to find both negative and positive spikes
    for k = 1:2
        if k == 2
            kdata = -hpdata; % flip the sign of the data to find positive spikes
        elseif k == 1
            kdata = hpdata;
        end
        
        % find peaks (spp) and troughs (spv) in the data
        [spp,spv] = FindPeaks(kdata);


        idx      = find(diff(spp) <= spkdur(2));       % find peak-to-peak durations within allowable range
        startdx  = spp(idx);
        startdx1 = spp(idx+1);

        % Loop over peaks
        for i = 1:length(startdx)
            spkmintic = spv((spv > startdx(i) & spv < startdx1(i))); % find the valley that is between the two peaks

            % If the height from valley to either peak is big enough, it could
            % be a spike
            max_height = max(abs(kdata(startdx1(i)) - kdata(spkmintic)),abs(kdata(startdx(i)) - kdata(spkmintic)));
            if max_height > thresh   % see if the peaks are big enough
                
                spikes(end+1,1) = spkmintic; % add timestamp to the spike list
                spikes(end,2)   = (startdx1(i)-startdx(i)); % add spike duration to list
                spikes(end,3)   = max_height;  % add spike amplitude to list

            end

        end
    end


    toosmall = [];
    toosharp = [];
    toobig = [];

    % now have all the info we need to decide if this thing is a spike or
    % not. Loop over spikes and subject criteria.
    for i = 1:size(spikes, 1)  % for each spike
        
        % re-define baseline to be period surrounding spike
        istart = max(1,round(spikes(i,1)-sur_time*fs));
        iend = min(length(hpdata),round(spikes(i,1)+sur_time*fs));
          
        alt_thresh = median(abs(hpdata(istart:iend)))*tmul;
        
        if spikes(i,3) > alt_thresh && spikes(i,3) > absthresh  % both parts together are bigger than thresh: so have some flexibility in relative sizes
            if spikes(i,2)*1000/fs > spkdur(1)    % spike wave cannot be too sharp: then it is either too small or noise
                if spikes(i,3) < too_high_abs
                    out(end+1,:) = spikes(i,:);  % add info of spike to output list
                    
                else
                    toobig(end+1) = spikes(i,1);
                end  
                
            else
                toosharp(end+1) = spikes(i,1);
            end
        else
            toosmall(end+1) = spikes(i,1);
        end
        
        
    end


    if ~isempty(out)
       
         % Re-align spikes to peak of the spikey component
         timeToPeak = [-.15,.15]; %Only look 150 ms before and 150 ms after the currently defined peak
         fullSurround = [-sur_time,sur_time]*fs;
         idxToPeak = timeToPeak*fs;
         
         for i = 1:size(out,1)
            currIdx = out(i,1);
            surround_idx = max(1,round(currIdx+fullSurround(1))):...
                min(round(currIdx+fullSurround(2)),length(hpdata));
            idxToLook = max(1,round(currIdx+idxToPeak(1))):...
                    min(round(currIdx+idxToPeak(2)),length(hpdata));  
            snapshot = data(idxToLook)-median(data(surround_idx)); % Look at the high frequency data (where the mean is substracted already)
            [~,I] = max(abs(snapshot)); % The peak is the maximum absolute value of this
            out(i,1) = idxToLook(1) + I - 1;
         end
        %}
    end
    %}



   all_spikes = [all_spikes;repmat(j,size(out,1),1) out];

   
end



gdf = all_spikes;
gdf = unique(gdf,'stable','rows');

%% sort by times and put ch first
if isempty(gdf) == 0
    gdf = sortrows(gdf,2); % sort by time

    %{
    times = gdf(:,1);
    chs = gdf(:,2);
    [times,I] = sort(times);
    chs = chs(I);
    gdf = [chs,times];
    %}
end



%% Remove those at beginning and end
if ~isempty(gdf)
    close_idx = close_to_edge*fs;
    gdf(gdf(:,2) < close_idx,:) = [];
    gdf(gdf(:,2) > size(values,1) - close_idx,:) = [];
end

%% remove duplicates
if ~isempty(gdf)
    keep = ones(size(gdf,1),1);

    % take diff of times
    diff_times = [inf;diff(gdf(:,2))];

    % take diff of chs
    diff_chs = [inf;diff(gdf(:,1))];

    % find those that are close in time and the same ch
    too_close = abs(diff_times) < 100e-3*fs & diff_chs == 0;

    keep(too_close) = 0;
    keep = logical(keep);

    n_removed = sum(~keep);
    gdf(~keep,:) = [];
end

%% Multichannel requirements
% Require spike to be on at least 2 channels and no more than half of the
% channels within 100 ms
if ~isempty(gdf)
    gdf =  multi_channel_requirements2(gdf,nchs,fs);
end



end