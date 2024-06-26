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
    options.tmul = 19; % minimum relative amplitude (compared to baseline). 
    options.absthresh = 100; % min absolute amplitude (uV). Can be a scalar or vector with one value per channel
    options.sur_time = 0.5; % surround time (in s) against which to compare for relative amplitude
    options.close_to_edge = 0.05; % time (in s) surrounding start and end of sample to ignore 
    options.too_high_abs = 1e3; % amplitude above which I reject it as artifact
    options.spkdur = [15 200]; % spike duration must be within this range (in ms)
    options.lpf1 = 30; % scalar: low pass cutoff freq for spikey component, or filter coefficient struct: lpf1.B, lpf1.A
    options.hpf  = 7; % scalar: high pass cutoff freq for spikey component, or filter cofficient struct: hpf.B, hpf.A
    options.multiChannelRequirements = true
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

% re-adjust the mean of the data to be zero
data = values - mean(values,1,'omitnan');

% Low pass filter to remove artifact
if ~isstruct(lpf1)
    all_nan_chs = all(isnan(data),1);
    temp_data = data(:,~all_nan_chs);
    lpdata = nan(size(data));
    lpdata_temp = eegfilt(temp_data, lpf1, 'lp',fs); % low pass filter
    lpdata(:,~all_nan_chs) = lpdata_temp;
else
    lpdata = filtfilt(lpf1.B,lpf1.A,data); % filter with passed inputs
end

% high pass filter to get the spikey part
if ~isstruct(hpf)
    all_nan_chs = all(isnan(lpdata),1);
    hpdata_all = nan(size(data));
    temp_data = lpdata(:,~all_nan_chs);
    hpdata_all_temp   = eegfilt(temp_data, hpf, 'hp',fs); % high pass filter
    hpdata_all(:,~all_nan_chs) = hpdata_all_temp;
else
    hpdata_all = filtfilt(hpf.B,hpf.A,lpdata);
end
    
% establish the baseline for the relative amplitude threshold
lthresh = median(abs(hpdata_all), 1); 
thresh  = lthresh*tmul;     % this is the final threshold we want to impose

% set absthresh for each channel
if length(absthresh) == 1, absthresh = repelem(absthresh, nchs); 
else, assert(length(absthresh) == nchs, 'absThresh must be a scalar or vector equal to # channels')
end

%% Iterate channels and detect spikes
for j = 1:nchs
    
    % initialize out array with final spike info
    out = [];
    
    hpdata = hpdata_all(:,j);

    % Skip if all nans
    if sum(isnan(hpdata)) > 0, continue; end
            
    % initialize array with tentative spike info
    spikes   = [];

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
        strtMean = mean([startdx, startdx1],2);

        % find the valley that is between the two peaks (get the index of
        % the valley that is closest to the mean of the peaks).
        [~,imn] = pdist2(spv, strtMean, 'euclidean', 'smallest', 1);
        spkmintic = spv(imn);

        % Loop over peaks
        for i = 1:length(startdx)
            % If the height from valley to either peak is big enough, it could
            % be a spike
            max_height = max(abs(kdata(startdx1(i)) - kdata(spkmintic(i))),abs(kdata(startdx(i)) - kdata(spkmintic(i))));
            if max_height > thresh(j)   % see if the peaks are big enough
                
                spikes(end+1,1) = spkmintic(i); % add timestamp to the spike list
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
        
        if spikes(i,3) > alt_thresh && spikes(i,3) > absthresh(j)  % both parts together are bigger than thresh: so have some flexibility in relative sizes
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
if ~isempty(gdf) && options.multiChannelRequirements
    gdf =  multi_channel_requirements2(gdf,nchs,fs);
end



end