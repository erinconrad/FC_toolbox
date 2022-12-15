function out = filter_canonical_freqs(values,fs)

%% Parameters
order = 4;
freqs = [1 4;...
    4 8;...
    8 12;...
    12 30;...
    30 70];

nfreqs = size(freqs,1);
out = nan(size(values,1),size(values,2),nfreqs);
for f = 1:nfreqs
    out(:,:,f) = bandpass_any(values,fs,freqs(f,:),order);
end

if 0
    figure
    ch = 7;
    set(gcf,'position',[10 10 1400 400])
    t = tiledlayout(1,nfreqs+1,'TileSpacing','tight','padding','tight');
    nexttile 
    plot(values(:,ch))
    for f = 1:nfreqs
        nexttile
        plot(out(:,ch,f))
        title(sprintf('%d-%d Hz',freqs(f,1),freqs(f,2)))
        set(gca,'fontsize',15)
    end
    title(t,sprintf('Bandpass IIR order %d filtfilt',order))
    set(gca,'fontsize',15)
end

end