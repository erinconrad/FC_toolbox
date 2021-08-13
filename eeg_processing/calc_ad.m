function ad_rat = calc_ad(values,fs)

nchs = size(values,2);
ad_rat = zeros(nchs,1);

for ich = 1:nchs
    X = values(:,ich);
    
    % turn any nans into the mean of the signal
    X(isnan(X)) = nanmean(X);
    
    % subtract mean
    X = X - mean(X);
    
    % Calculate fft
    Y = fft(X);

    % Get power
    P = abs(Y).^2;
    freqs = linspace(0,fs,length(P)+1);
    freqs = freqs(1:end-1);

    % Take first half
    P = P(1:ceil(length(P)/2));
    freqs = freqs(1:ceil(length(freqs)/2));

    %plot(freqs,P);

    alpha = sum(P(freqs>=8 & freqs<=13));
    delta = sum(P(freqs>=1 & freqs<=4));
    
    ad_rat(ich) = alpha/delta;
    
    if 0
        figure
        tiledlayout(1,2)
        nexttile
        plot(X)
        
        nexttile
        plot(freqs,P)
        hold on
        plot([1 1],ylim,'b')
        plot([4 4],ylim,'b')
        plot([8 8],ylim,'r')
        plot([13 13],ylim,'r')
        title(sprintf('Ratio = %1.2f',alpha/delta));
        pause
        close(gcf)
        
    end
    
end


end