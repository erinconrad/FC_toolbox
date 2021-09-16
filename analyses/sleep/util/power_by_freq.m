function [P,freqs] = power_by_freq(thing,fs)

X = thing - nanmean(thing);
X(isnan(X)) = nanmean(X);
Y = fft(X);
P = (abs(Y)).^2; % power
freqs = linspace(0,fs,length(P)+1);
freqs = freqs(1:end-1);
% Take first half
P = P(1:ceil(length(P)/2));
freqs = freqs(1:ceil(length(freqs)/2));

end