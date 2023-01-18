function se_all = spectral_entropy(values,fs,tw,do_tw)

nchs = size(values,2);

if do_tw
    % divide into time windows
    tw = round(tw*fs);
    times = 1:tw:size(values,1);
    
    se_all = nan(nchs,length(times)-1);
    
    for t = 1:length(times)-1
        for ich = 1:nchs
            xn = values(times(t):times(t+1),ich);
            xn(isnan(xn)) = nanmean(xn);
            if all(isnan(xn))
                continue
            end
            [se,te] = pentropy(xn,fs);
        
            if 0
                tt = 1/fs*(0:length(xn)-1);
                figure
                nexttile
                plot(tt,xn)
                nexttile
                plot(te,se)
        
            end
            se_all(ich,t) = nanmean(se);
        
        
        end
    end
    
    se_all = nanmean(se_all,2);

else

    se_all = nan(nchs,1);
    
    for ich = 1:nchs
        xn = values(:,ich);
        xn(isnan(xn)) = nanmean(xn);
        if all(isnan(xn))
            continue
        end
        [se,te] = pentropy(xn,fs);
    
        if 0
            t = 1/fs*(0:length(xn)-1);
            figure
            nexttile
            plot(t,xn)
            nexttile
            plot(te,se)
    
        end
        se_all(ich) = nanmean(se);
    
    
    end

end




end