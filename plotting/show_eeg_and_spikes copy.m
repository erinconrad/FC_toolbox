function show_eeg_and_spikes(values,chLabels,gdf,fs)

figure
set(gcf,'position',[62 104 1145 701])

offset = 0;
ch_offsets = zeros(size(values,2),1);
ch_bl = zeros(size(values,2),1);
dur = size(values,1)/fs;

for ich = 1:size(values,2)

    plot(linspace(0,dur,size(values,1)),values(:,ich)-offset,'k');
    
    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(values(:,ich));
    hold on
    text(dur+0.05,ch_bl(ich),sprintf('%s',chLabels{ich}))
    
    if ich<size(values,2)
        if ~isnan(min(values(:,ich)) - max(values(:,ich+1)))
            offset = offset - (min(values(:,ich)) - max(values(:,ich+1)));
        end
    end
end

for s = 1:size(gdf,1)
    %index = spikes(s,1);
    index = gdf(s,2);
    
    % convert index to time
    time = index/fs;
    
    ch = gdf(s,1);
    offset_sp = ch_offsets(ch);
    
    value_sp = values(round(index),ch);
    
    plot(time,value_sp - offset_sp,'ro')
    
end


pause
close(gcf)

end