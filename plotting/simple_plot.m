function simple_plot(tout,out,chs,im)


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
    chs = find(is_run);
end

offset = 0;
ch_offsets = zeros(length(chs),1);
ch_bl = zeros(length(chs),1);
for i = 1:length(chs)
    ich = chs(i);
    plot(linspace(0,dur,size(values,1)),values(:,ich)-offset,'linewidth',2);
    hold on
    ch_offsets(i) = offset;
    ch_bl(i) = -offset + nanmedian(values(:,ich));

    text(dur+0.05,ch_bl(i),sprintf('%s',labels{ich}),'fontsize',20)

    if i<length(chs)
        if ~isnan(min(values(:,ich)) - max(values(:,chs(i+1))))
            offset = offset - (min(values(:,ich)) - max(values(:,chs(i+1))));
        end
    end
end
xlabel('Time (seconds)')
set(gca,'fontsize',20)
    
    



end