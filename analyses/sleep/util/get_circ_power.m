function P_in_bin = get_circ_power(rate,times)

pbin = [22 26];

% Get frequency in Hz of the 24ish hour band I want
pbin = pbin*3600;
fbin = 1./pbin;

fs = 1/median(diff(times)); % how often I am spacing this in seconds
X = rate - nanmean(rate); % subtract dc component
X(isnan(X)) = nanmean(X); % set nans to mean
Y = fft(X); % fft
P = (abs(Y)).^2; % power
freqs = linspace(0,fs,length(P)+1);
freqs = freqs(1:end-1);
% Take first half
P = P(1:ceil(length(P)/2));
freqs = freqs(1:ceil(length(freqs)/2));
Pnorm = P;

% find freqs in the fbin
freqs_in_bin = freqs > fbin(2) & freqs < fbin(1);

% Get summed relative power in that bin
P_in_bin = sum(Pnorm(freqs_in_bin))/sum(P);

%period = 1./freqs;

end