function all_plv = plv_calc(values,fs)

%% Define frequency bands
freqs = [1 4;...
    4 8;...
    8 12;...
    12 30;...
    30 70
    1 70];


nfreqs = size(freqs,1);
nchs = size(values,2);


%% initialize output vector
all_plv = nan(nchs,nchs,nfreqs);

for ich = 1:nchs
    curr_values = values(:,ich);
    curr_values(isnan(curr_values)) = nanmean(curr_values);
    values(:,ich) = curr_values;
end

%% Remove nan rows, keeping track of which I am removing
nan_rows = any(isnan(values),1); % find channels with nans for any time points
values_no_nans = values(:,~nan_rows);
nchs_no_nans = size(values_no_nans,2);
temp_plv = nan(nchs_no_nans,nchs_no_nans,nfreqs);
A = values_no_nans;
nchs = size(A,2);

for f = 1:nfreqs
    filtSpec.range = freqs(f,:);
    % Want 4-5 cycles. ncycles = T x freq_ramge = order/fs * freq-range
    % So order should be about ncycles*fs/freq_range
    filtSpec.order = round(fs*5/(min(filtSpec.range)));
    filtPts = fir1(filtSpec.order, 2/fs*filtSpec.range);
    filteredData = filter(filtPts, 1, A, [], 2);

    % Get phase of each signal
    phase = nan(size(A));
    for ich = 1:nchs
        phase(:,ich)= angle(hilbert(filteredData(:,ich)));
    end

    % Get PLV
    plv = nan(nchs,nchs);
    for ich = 1:nchs
        for jch = ich+1:nchs
            e = exp(1i*(phase(:,ich) - phase(:,jch)));
            plv(ich,jch) = abs(sum(e,1))/size(phase,1);
            plv(jch,ich) = abs(sum(e,1))/size(phase,1);
        end
    end
    temp_plv(:,:,f) = plv;

    if 0
        figure
        nexttile
        plot(filteredData(:,2))
        hold on
        plot(filteredData(:,3))
        nexttile
        plot(phase(:,2))
        hold on
        plot(phase(:,3))
    end
end

%% Put the non-nans back
all_plv(~nan_rows,~nan_rows,:) = temp_plv;
all_plv(logical(repmat(eye(nchs,nchs),1,1,nfreqs))) = nan;

if 0
    figure; tiledlayout(1,6)
    for i = 1:6
        nexttile
        turn_nans_gray(all_plv(:,:,i))
        colorbar
    end
end

end