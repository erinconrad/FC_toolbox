function all_power = bp_calc(run_values,fs,skip)

freqs = [0.5 4;...
    4 8;...
    8 12;...
    12 30;...
    30 80];
nfreqs = size(freqs,1);

nchs = size(run_values,2);
all_power = nan(nchs,nfreqs);

% turn skip channels into zeros
run_values(:,skip) = 0;

for f = 1:nfreqs
    curr_freq = freqs(f,:);
    p = bandpower(run_values,fs,curr_freq);
    all_power(:,f) = p;

end

all_power(skip,:) = nan;


end