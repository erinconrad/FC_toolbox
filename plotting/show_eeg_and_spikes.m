function show_eeg_and_spikes(values,chLabels,gdf,fs,which_chs)

%figure
%set(gcf,'position',[62 104 1145 701])

offset = 0;
%ch_offsets = zeros(size(values,2),1);
%ch_bl = zeros(size(values,2),1);
ch_offsets = zeros(length(which_chs),1);
ch_bl = zeros(length(which_chs),1);
dur = size(values,1)/fs;
%dur = size(values,1);

for i = 1:length(which_chs)
    ich = which_chs(i);
    plot(linspace(0,dur,size(values,1)),values(:,ich)-offset,'k','linewidth',2);
    %plot(values(:,ich)-offset,'k','linewidth',2);
    
    ch_offsets(i) = offset;
    ch_bl(i) = -offset + nanmedian(values(:,ich));
    hold on
    text(dur+0.05,ch_bl(i),sprintf('%s',chLabels{ich}))
    
    if i<length(which_chs)
        if ~isnan(min(values(:,ich)) - max(values(:,which_chs(i+1))))
            offset = offset - (min(values(:,ich)) - max(values(:,which_chs(i+1))));
        end
    end
end

for s = 1:size(gdf,1)
    %index = spikes(s,1);
    index = gdf(s,2);
    
    % convert index to time
    time = index/fs;
    %time = index;
    
    ch = gdf(s,1);
    if ~ismember(ch,which_chs)
        continue
    end
    offset_sp = ch_offsets(ismember(which_chs,ch));
    
    value_sp = values(round(index),ch);
    
    plot(time,value_sp - offset_sp,'bo','linewidth',2)
    
end

grid off
xticklabels([])
axis off
%pause
%close(gcf)

end