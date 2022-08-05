function all_coherence = faster_coherence_calc(values,fs)

%% Parameters
window = fs * 2;
NFFT = window;
freqs = [1 4;...
    4 8;...
    8 12;...
    12 30;...
    30 70
    1 70];
nfreqs = size(freqs,1);

nchs = size(values,2);

%% initialize output vector
all_coherence = nan(nchs,nchs,nfreqs);

for ich = 1:nchs
    curr_values = values(:,ich);
    curr_values(isnan(curr_values)) = nanmean(curr_values);
    values(:,ich) = curr_values;
end


%% Remove nan rows, keeping track
nan_rows = any(isnan(values),1); % find channels with nans for any time points
values_no_nans = values(:,~nan_rows);
nchs_no_nans = size(values_no_nans,2);
temp_coherence = nan(nchs_no_nans,nchs_no_nans,nfreqs);


for ich = 1:nchs_no_nans

    % Do MS cohere on full thing
    [cxy,f] = mscohere(values_no_nans(:,ich),values_no_nans,hamming(window),[],NFFT,fs);
    %[cxy,f] = mscohere(values(:,ich),values,hamming(window),[],NFFT,fs);

    % Average coherence in frequency bins of interest
    for i_f = 1:nfreqs
        temp_coherence(:,ich,i_f) = ...
            nanmean(cxy(f >= freqs(i_f,1) & f <= freqs(i_f,2),:),1);

        temp_coherence(:,ich,i_f) = ...
            nanmean(cxy(f >= freqs(i_f,1) & f <= freqs(i_f,2),:),1);
    end
    
end

%% Put the non-nans back
all_coherence(~nan_rows,~nan_rows,:) = temp_coherence;
all_coherence(logical(repmat(eye(nchs,nchs),1,1,nfreqs))) = nan;


end