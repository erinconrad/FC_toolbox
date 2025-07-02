function gdf = clean_detector_one_channel_fast(data,fs)
    
    %{
    Inputs:
    - data: an nsamples x 1 array of eeg data (only operates on one channel).
    - fs: sampling rate in Hz
    
    Outputs:
    - gdf: an nspike x 5 array with data for each detected spike. Spikes are in
    ascending order of time and columns are the following
       - column 1: spike channel
       - column 2: spike time in SAMPLES (not s)
       - column 3: spike duration in SAMPLES
       - column 4: spike amplitude in uV
       - column 5: spike sign (1=positive spike, -1=negative spike)
    
    Dependencies:
    - eegfilt
    - FindPeaks
    
    Information:
    - This was originally written by Camilo Bermudez 7/31/13 and modified by
    Erin Conrad at UPenn 12/9/22.
    - Modified by M.A. Kramer 07/2025 to speed up.
        x. Only operates on 1 channel.
    %}
    
    %%%% Parameters
    params.tmul            = 19;       % minimum relative amplitude (compared to baseline)
    params.absthresh       = 100;      % minimum absolute amplitude (uV)
    params.sur_time        = 0.5;      % surround time (in s) against which to compare for relative amplitude
    params.close_to_edge   = 0.05;     % time (in s) surrounding start and end of sample to ignore
    params.too_high_abs    = 1e3;      % amplitude above which I reject it as artifact
    params.spkdur          = [15 200]; % spike duration must be within this range (in ms)
    params.spkdur          = params.spkdur*fs/1000;   % convert above to samples;
    params.lpf1            = 30;       % low pass filter for spikey component
    params.hpf             = 7;        % high pass filter for spikey component

    %%%% Initialize things
    all_spikes  = [];

    % initialize out array with final spike info
    out = [];

    % Skip if all nans
    %if sum(isnan(data)) > 0, continue; end

    % re-adjust the mean of the data to be zero
    data = data - nanmean(data);

    % Low pass filter to remove artifact
    lpdata = eegfilt(data, params.lpf1, 'lp',fs); % low pass filter

    % high pass filter to get the spikey part
    hpdata = eegfilt(lpdata, params.hpf, 'hp',fs); % high pass filter
    
    % establish the baseline for the relative amplitude threshold
    lthresh = median(abs(hpdata));
    thresh  = lthresh*params.tmul;     % this is the final threshold we want to impose
    
    %%%% Run the spike detector to find both negative and positive spikes %

    %%%% Positive spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kdata = hpdata;

    % find peaks (spp) and troughs (spv) in the data
    [spp,spv] = FindPeaks(kdata);
    idx       = find(diff(spp) <= params.spkdur(2));                       % find peak-to-peak durations within allowable range
    startdx   = spp(idx);
    startdx1  = spp(idx+1);

    % Loop over peaks
    spikes   = zeros(length(startdx),4);                                   % initialize array with spike info
    counter  = 1;
    for i=1:length(startdx)
        ktemp      = kdata(startdx(i):startdx1(i));                        % Get the data between the two peaks.
        max_height = max([ktemp(1)-min(ktemp), ktemp(end)-min(ktemp)]);    % Compute distance between each peak and min between them.
        if max_height > thresh                                             % see if the peaks are big enough
            [~,min_index]       = min(ktemp);                              % ... if so, store it.
            spkmintic           = startdx(i) + min_index;
            spikes(counter,1)   = spkmintic;                               % add timestamp to the spike list
            spikes(counter,2)   = (startdx1(i)-startdx(i));                % add spike duration to list
            spikes(counter,3)   = max_height;                              % add spike amplitude to list
            spikes(counter,4)   = 1;                                       % add spike sign to list
            counter=counter+1;
        end
    end
    positive_spikes = spikes(1:counter-1,:);

    %%%% Negative spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kdata = -hpdata;

    % find peaks (spp) and troughs (spv) in the data
    [spp,spv] = FindPeaks(kdata);
    idx       = find(diff(spp) <= params.spkdur(2));       % find peak-to-peak durations within allowable range
    startdx   = spp(idx);
    startdx1  = spp(idx+1);

    % Loop over peaks
    % initialize array with tentative spike info
    spikes   = zeros(length(startdx),3);
    counter  = 1;
    for i=1:length(startdx)
        ktemp      = kdata(startdx(i):startdx1(i));                        % Get the data between the two peaks.
        max_height = max([ktemp(1)-min(ktemp), ktemp(end)-min(ktemp)]);    % Compute distnace between each peak and min between them.
        if max_height > thresh                                             % see if the peaks are big enough
            [~,min_index]       = min(ktemp);
            spkmintic           = startdx(i) + min_index;
            spikes(counter,1)   = spkmintic;                               % add timestamp to the spike list
            spikes(counter,2)   = (startdx1(i)-startdx(i));                % add spike duration to list
            spikes(counter,3)   = max_height;                              % add spike amplitude to list
            spikes(counter,4)   = -1;                                      % add spike sign to list
            counter=counter+1;
        end
    end
    negative_spikes = spikes(1:counter-1,:);

    %%%% Combine the positive and negative spikes. %%%%%%%%%%%%%%%%%%%%%%%%
    spikes   = [positive_spikes; negative_spikes];

    toosmall = [];
    toosharp = [];
    toobig   = [];

    % now have all the info we need to decide if this thing is a spike or
    % not. Loop over spikes and subject criteria.
    for i = 1:size(spikes, 1)  % for each spike

        % re-define baseline to be period surrounding spike
        istart = max(1,round(spikes(i,1)-params.sur_time*fs));
        iend = min(length(hpdata),round(spikes(i,1)+params.sur_time*fs));

        alt_thresh = median(abs(hpdata(istart:iend)))*params.tmul;

        if spikes(i,3) > alt_thresh && spikes(i,3) > params.absthresh  % both parts together are bigger than thresh: so have some flexibility in relative sizes
            if spikes(i,2)*1000/fs > params.spkdur(1)    % spike wave cannot be too sharp: then it is either too small or noise
                if spikes(i,3) < params.too_high_abs
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
        fullSurround = [-params.sur_time,params.sur_time]*fs;
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

    end

    all_spikes = [all_spikes;repmat(1,size(out,1),1) out];

    gdf = all_spikes;
    gdf = unique(gdf,'stable','rows');

    %%%% sort by times and put ch first
    if isempty(gdf) == 0
        gdf = sortrows(gdf,2); % sort by time
    end

    %%%% Remove those at beginning and end
    if ~isempty(gdf)
        close_idx = params.close_to_edge*fs;
        gdf(gdf(:,2) < close_idx,:) = [];
        gdf(gdf(:,2) > size(data,1) - close_idx,:) = [];
    end

    %%%% remove duplicates
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

end