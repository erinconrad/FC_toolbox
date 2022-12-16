function re = relative_entropy(values,fs)

%% Get filtered signal
out = filter_canonical_freqs(values,fs);
nfreqs = size(out,3);
nchs = size(values,2);

re = nan(nchs,nchs,nfreqs);


% loop over frequences
for f = 1:nfreqs

    filteredData = out(:,:,f);

    for ich = 1:nchs
        for jch = ich+1:nchs
            h1 = histcounts(filteredData(:,ich),10);
            h2 = histcounts(filteredData(:,jch),10);
            h1 = h1/sum(h1); h2 = h2/sum(h2); % normalize?
            S1 = sum(h1*log(h1/h2));
            S2 = sum(h2*log(h2/h1));
            re(ich,jch,f) = max([S1,S2]);
            re(jch,ich,f) = re(ich,jch,f);
        end
    end
end

if 0
    figure
    for i = 1:nfreqs
    nexttile
    turn_nans_gray(re(:,:,i))
    end

end


end