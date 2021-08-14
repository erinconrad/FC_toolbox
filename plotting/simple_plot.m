function simple_plot(tout,out,chs,im,gdf,only_run,bad)


fs = out.fs;
clean_labels = out.clean_labels;

if ~isempty(chs)
    new_chs = nan(length(chs),1);
        
    % get ch index
    for j = 1:length(chs)
        curr_ch = chs{j};

        % find index
        curr_idx = find(strcmp(curr_ch,clean_labels));

        new_chs(j) = curr_idx;
    end

    chs = new_chs;
end


figure
set(gcf,'position',[10 10 1200 800])
  
values = tout.montage(im).values;
labels = out.montage(im).labels;
is_run = out.montage(im).is_run;
dur = size(values,1)/fs;

if isempty(chs)
    if only_run
        chs = find(is_run);
        
        % Convert spike channels
        all_chs = 1:length(labels);
        
    else
    
        chs = 1:length(labels);
    end
end

offset = 0;
ch_offsets = zeros(length(chs),1);
ch_bl = zeros(length(chs),1);
last_min = nan;
for i = 1:length(chs)
    ich = chs(i);
    
    if sum(~isnan(values(:,ich))) ~=0
        
        if ismember(ich,bad)
            plot(linspace(0,dur,size(values,1)),values(:,ich)-offset,'r');
        else
            plot(linspace(0,dur,size(values,1)),values(:,ich)-offset,'k');
        end
        hold on
        ch_offsets(i) = offset;
        ch_bl(i) = -offset + nanmedian(values(:,ich));

        text(dur+0.05,ch_bl(i),sprintf('%s',labels{ich}),'fontsize',20)
        
        last_min = min(values(:,ich));
    end

   
    if i<length(chs)
        
        
        if ~isnan(max(values(:,chs(i+1)))) & ~isnan(last_min)
            offset = offset - (last_min - max(values(:,chs(i+1))));
        end
            %{
        if ~isnan(min(values(:,ich)) - max(values(:,chs(i+1))))
            offset = offset - (min(values(:,ich)) - max(values(:,chs(i+1))));
        end
            %}
    end
    
end
xlabel('Time (seconds)')
set(gca,'fontsize',20)
    
for s = 1:size(gdf,1)
    %index = spikes(s,1);
    index = gdf(s,2);
    
    % convert index to time
    time = index/fs;
    
    ch = gdf(s,1);
    offset_sp = ch_offsets(ch);
    
    value_sp = values(round(index),ch);
    
    plot(time,value_sp - offset_sp,'ko','markersize',10,'linewidth',2)
    
end



end