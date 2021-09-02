function gdf = detector_alt(values,fs)

%{
This is the spike detection algorithm used for the implant effect paper.
%}

tmul = 19;
absthresh = 100;
sur_time = 0.5;
close_to_edge = 0.05;

% Initialize parameters
too_high_abs = 1e3; % tmul above which I reject it as artifact
spkdur = [15 200];                % spike duration must be less than this in ms
spkdur = spkdur*fs/1000;   % convert to points;
lpf1     = 30; % low pass filter for spikey component
hpf  = 7; % high pass filter for spikey component
vlpf = 0.5; % very low pass filter to exclude spikes if very low freq shift

% Initialize things
all_spikes  = [];
nchs = size(values,2);

% Iterate channels and detect spikes
for j = 1:nchs
    
    out = [];
    
    % specify ch
    dd = j;
    data = values(:,dd);
    
    if sum(isnan(data)) > 0, continue; end
    
    % re-adjust the mean of the data to be zero
    data = data - nanmean(data);
        
    spikes   = [];

    % Low pass filter to remove artifact
    lpdata = eegfilt(data, lpf1, 'lp',fs); % low pass filter
    % first look at the high frequency data for the 'spike' component
    hpdata   = eegfilt(lpdata, hpf, 'hp',fs); % high pass filter
    %hf_values(:,dd) = HFdata;

    if 0
        figure
        tiledlayout(3,1)
        nexttile
        plot(data)
        nexttile
        plot(lpdata)
        nexttile
        plot(hpdata)
    end

    lthresh = median(abs(hpdata)); 
    thresh  = lthresh*tmul;     % this is the final threshold we want to impose

    for k = 1:2
        if k == 2
            %kdata = -lpdata;
            kdata = -hpdata;
        else
            %kdata = lpdata;
            kdata = hpdata;
        end
        
        [spp,spv] = FindPeaks(kdata);

        idx      = find(diff(spp) <= spkdur(2));       % find the durations less than or equal to that of a spike
        startdx  = spp(idx);
        startdx1 = spp(idx+1);

        % Loop over peaks
        for i = 1:length(startdx)
            spkmintic = spv((spv > startdx(i) & spv < startdx1(i))); % find the valley that is between the two peaks

            % If the height from valley to either peak is big enough, it could
            % be a spike
            max_height = max(abs(kdata(startdx1(i)) - kdata(spkmintic)),abs(kdata(startdx(i)) - kdata(spkmintic)));
            if max_height > thresh   % see if the peaks are big enough
                
                spikes(end+1,1) = spkmintic;                                  % add timestamp to the spike list
                spikes(end,2)   = (startdx1(i)-startdx(i))*1000/fs;         % add spike duration to list
                spikes(end,3)   = max_height;    % add spike amplitude to list

            end

        end
    end


  
    if 0
        figure
        plot(data)
        hold on
        plot(spikes(:,1),data(spikes(:,1)),'o')

    end

    tooshort = [];
    toosmall = [];
    toosharp = [];
    toobig = [];

    % now have all the info we need to decide if this thing is a spike or not.
    for i = 1:size(spikes, 1)  % for each spike
        
        % re-define baseline to be 2 seconds surrounding
        surround = sur_time;
        istart = max(1,round(spikes(i,1)-surround*fs));
        iend = min(length(hpdata),round(spikes(i,1)+surround*fs));
          
        alt_thresh = median(abs(hpdata(istart:iend)))*tmul;
        
        if spikes(i,3) > alt_thresh && spikes(i,3) > absthresh            % both parts together are bigger than thresh: so have some flexibility in relative sizes
            if spikes(i,2) > spkdur(1)     % spike wave cannot be too sharp: then it is either too small or noise
                if spikes(i,3) < too_high_abs
                    out(end+1,1) = spikes(i,1);         % add timestamp of spike to output list
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
       
        %get_spike_details(out,data,fndata,HFdata,fs)
        %
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
            snapshot = hpdata(idxToLook)-median(hpdata(surround_idx)); % Look at the high frequency data (where the mean is substracted already)
            [~,I] = max(abs(snapshot)); % The peak is the maximum absolute value of this
            out(i,1) = idxToLook(1) + I - 1;
         end
        %}
    end
    %}



   all_spikes = [all_spikes;out,repmat(dd,length(out),1)];

    if 0
        figure
        tiledlayout(3,1)
        nexttile
        plot(data)
        hold on
        if ~isempty(out)
            plot(out(:,1),data(out(:,1)),'o')
        end
        nexttile
        plot(lpdata)
        nexttile
        plot(hpdata)
        

    end
   
end



gdf = all_spikes;
gdf = unique(gdf,'stable','rows');

if isempty(gdf) == 0
    times = gdf(:,1);
    chs = gdf(:,2);
    [times,I] = sort(times);
    chs = chs(I);
    gdf = [chs,times];
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

% Multichannel requirements
if ~isempty(gdf)
    gdf =  multi_channel_requirements(gdf,nchs,fs);
end



end