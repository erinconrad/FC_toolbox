function show_eeg_and_spikes_select(values,gdf,fs,to_plot,show_spikes)
color = [0.8500, 0.3250, 0.0980];
%figure
%set(gcf,'position',[62 104 1145 701])

nplot = sum(to_plot);
offset = 0;
ch_offsets = zeros(nplot,1);
ch_idx = nan(nplot,1);

count = 0;
for ich = 1:size(values,2)
    
    if to_plot(ich) == 1
        count = count+1;

        plot(values(:,ich)-offset,'k','linewidth',2);
    
        ch_offsets(count) = offset;
        ch_idx(count) = ich;
        hold on
    %text(dur+0.05,ch_bl(ich),sprintf('%s',chLabels{ich}))
    
        % Add spikes
        if show_spikes
            [lia,loc] = ismember(ich,gdf(:,1));
            if lia
                for k = 1:length(loc)
                    sp = gdf(loc,:);
                    plot(sp(:,2),values(sp(:,2),ich)-offset,'o','color',color,...
                        'linewidth',3,'markersize',10);
                end
            end
        end
            
        if ich<size(values,2)
            for k = ich+1:size(values,2)
                if to_plot(k)
                    if ~isnan(min(values(:,ich)) - max(values(:,k)))
                        offset = offset - (min(values(:,ich)) - max(values(:,k)));
                    end
                    break
                end
                
            end
        end
    end
end



%pause
%close(gcf)

end