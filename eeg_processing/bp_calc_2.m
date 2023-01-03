function [all_power,rel_power] = bp_calc_2(run_values,fs,skip,tw,do_tw)

freqs = get_frequencies; 
freqs = freqs(2:end,:);
nfreqs = size(freqs,1);

nchs = size(run_values,2);
all_power = nan(nchs,nfreqs);

% turn skip channels into zeros
run_values(:,skip) = 0;
run_values(:,all(isnan(run_values),1)) = 0;

if do_tw

    % divide into time windows
    tw = round(tw*fs);
    times = 1:tw:size(run_values,1);
    
    for f = 1:nfreqs
        curr_freq = freqs(f,:);
        all_p = nan(length(times)-1,size(run_values,2));
        for t = 1:length(times)-1
            p = bandpower(run_values(times(t):times(t+1),:),fs,curr_freq);
            all_p(t,:) = p;
        end
        all_power(:,f) = nanmean(all_p,1);
    
    end

else
    for f = 1:nfreqs
        curr_freq = freqs(f,:);
        p = bandpower(run_values,fs,curr_freq);
        all_power(:,f) = p;
    
    end
end

all_power(skip,:) = nan;
all_power(all(isnan(run_values),1),:) = nan;


% also get relative power (divide signal by power in whole thing (first freq band is broadband)
rel_power = all_power./all_power(:,1); 


end