function re = relative_entropy(values,fs)

%% Define frequency bands
freqs = [1 4;...
    4 8;...
    8 12;...
    12 30;...
    30 70;...
    1 70];

nfreqs = size(freqs,1);
nchs = size(values,2);
re = nan(nchs,nchs,nfreqs);

%% Remove nans
for ich = 1:nchs
    curr_values = values(:,ich);
    curr_values(isnan(curr_values)) = nanmean(curr_values);
    values(:,ich) = curr_values;
end

% loop over frequences
for f = 1:nfreqs

    % filter
    filtSpec.range = freqs(f,:);
    % Want 4-5 cycles. ncycles = T x freq_ramge = order/fs * freq-range
    % So order should be about ncycles*fs/freq_range
    filtSpec.order = round(fs*5/(min(filtSpec.range)));
    filtPts = fir1(filtSpec.order, 2/fs*filtSpec.range);
    filteredData = filter(filtPts, 1, values, [], 2);


    for ich = 1:nchs
        for jch = ich+1:nchs
            h1 = histcounts(filteredData(:,ich),10);
            h2 = histcounts(filteredData(:,jch),10);
            h1 = h1/sum(h1); h2 = h2/sum(h2); % normalize?
            S1 = sum(h1*log(h1/h2));
            S2 = sum(h2*log(h2/h1));
            re(ich,jch,f) = max([S1,S2]);
            re(jch,ich,f) = re(ich,jch,f);
        end
    end
end

if 0
    figure
    nexttile
    plot(values(:,1))
    hold on
    plot(values(:,2))
    nexttile
    plot(h1)
    hold on
    plot(h2)
    nexttile
    turn_nans_gray(re(:,:,1))

end


end