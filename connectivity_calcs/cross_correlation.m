function [mb,lb] = cross_correlation(values,fs)


max_lag = 200; % 200 ms
ml = round(max_lag * 1e-3 *fs);

% Calculate cross correlation. 'Normalized' normalizes the sequence so that
% the autocorrelations at zero lag equal 1
[r,lags] = xcorr(values,ml,'normalized');

% find the maximum xcorr, and its corresponding lag, for each channel pair
[M,I] = max(r,[],1);
nlags = lags(I);

% back out what channels are what
mb = reshape(M,size(values,2),size(values,2));
lb = reshape(nlags,size(values,2),size(values,2));

% Make anything that is a nan in mb a nan in lb
lb(isnan(mb)) = nan;

% make anything at the edge of lags nan in lb
lb(lb==ml) = nan;
lb(lb==-ml) = nan;

lb = lb/fs;



if 0
    figure
    set(gcf,'position',[10 10 1400 400])
    nexttile
    turn_nans_gray(mb)
    nexttile
    turn_nans_gray(lb)
    nexttile
    turn_nans_gray(corr(values))
    
end


end