function show_scalp_eeg(values,fs,labels)

%% Parameters
break_points = {'F4-C4';'F4-F8'};
secs_per_plot = 15;

%% Initialize fig
dur = size(values,1)/fs;
nchs = length(labels);
nwindows = round(dur/secs_per_plot);
dur_indices = round(size(values,1)/nwindows);
if nwindows ~= 4, error('what'); end

figure
set(gcf,'position',[10 10 1200 800])
tiledlayout(nwindows,1,'tilespacing','tight')

for t = 1:nwindows
    
    nexttile
    
    tindices = (t-1)*dur_indices+1:t*dur_indices;
    eeg = values(tindices,:);
    
    offset = 0;
    ch_offsets = zeros(nchs,1);
    ch_bl = zeros(nchs,1);

    for ich = 1:nchs

        plot(linspace(0,secs_per_plot,size(eeg,1)),eeg(:,ich) - offset,'k');
        hold on
        ch_offsets(ich) = offset;
        ch_bl(ich) = -offset + nanmedian(eeg(:,ich));

        text(secs_per_plot+0.05,ch_bl(ich),sprintf('%s',labels{ich}),'fontsize',20)
        
        if ich < nchs

            if any(ismember(labels{ich},break_points))
                offset = offset - 2*(min(eeg(:,ich)) - max(eeg(:,ich+1)));
            else
                offset = offset - (min(eeg(:,ich)) - max(eeg(:,ich+1)));
            end

        end

    end
    
    if t == nwindows
        xlabel('Time (seconds)')
    else
        xticklabels([])
    end
    set(gca,'fontsize',15)
    
end



end