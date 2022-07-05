function all_coherence = faster_coherence_calc(values,fs)

%% Parameters
hws = 0.5; %2
overlaps = 0.3; %1
freqs = [0.5 4;...
    4 8;...
    8 12;...
    12 30;...
    30 80];
nfreqs = size(freqs,1);

nchs = size(values,2);

% convert hamming window and overlap from seconds to segments
hw = round(hws*fs);
overlap = round(overlaps*fs);

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

%% Second matrix

%% Do MS cohere on full thing
[cxy,f] = mscohere(values_no_nans,values_no_nans,hw,overlap,[],fs,'mimo');

%% Average coherence in frequency bins of interest
for i_f = 1:nfreqs
    temp_coherence(:,:,i_f) = ...
        nanmean(cxy(f >= freqs(i_f,1) & f <= freqs(i_f,2),:,:),1);

    temp_coherence(:,:,i_f) = ...
        nanmean(cxy(f >= freqs(i_f,1) & f <= freqs(i_f,2),:,:),1);
end


end