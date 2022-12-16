function se_all = spectral_entropy(values,fs)

nchs = size(values,2);
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