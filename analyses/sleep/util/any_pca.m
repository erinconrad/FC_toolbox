function any_pca(spike_bins,times,locs,names,title_text)

%ex_pts = [21,55,23,13,50,77];

npts = size(spike_bins,1);
nbins = size(spike_bins,2);
orig_bins = spike_bins;
% Subtract mean
spike_bins = (spike_bins - nanmean(spike_bins,2))./nanstd(spike_bins,[],2);

[coeff,score,latent] = pca(spike_bins,'Rows','pairwise');

% Top scorers
[~,top] = max(score,[],1);

% Bottom scorers
[~,bottom] = min(score,[],1);


top_12 = [top(1:2),bottom(1:2)];

figure
set(gcf,'position',[10 10 1300 600])
t = tiledlayout(2,6,'tilespacing','tight','padding','compact');

nexttile([1 2])
stem(latent,'linewidth',2)
title('Principal component variances')
set(gca,'fontsize',15)
xlabel('Component')
ylabel('Variance')

nexttile([1 2])
p1= plot(times,coeff(:,1),'linewidth',2);
hold on
p2 = plot(times,coeff(:,2),'linewidth',2);
plot([0 0],ylim,'k--','linewidth',2)
%p3 = plot(times,coeff(:,3),'linewidth',2);
ylabel('Coefficients')
title('Top principal components')
legend([p1 p2],{'Component 1','Component 2'},'fontsize',15,'location','northwest')
set(gca,'fontsize',15)
xlabel('Hours')
xlim([times(1) times(end)])

%{
nexttile
plot(score(:,1),score(:,2),'o','linewidth',2)
%scatter3(score(:,1),score(:,2),score(:,3))
xlabel('Score 1st component')
ylabel('Score 2nd component')
title('Scores across patients')
%zlabel('Score 3nd component')
set(gca,'fontsize',15)
%}

nexttile([1 2])
c = 1;
tloc = strcmp(locs,'temporal');
oloc = strcmp(locs,'other');
plot(1+0.05*randn(sum(tloc),1),score(tloc,c),'o','linewidth',2)
hold on
plot(2+0.05*randn(sum(oloc),1),score(oloc,c),'o','linewidth',2)
xlim([0.5 2.5])
plot(xlim,[0 0],'k--','linewidth',2)
p = ranksum(score(tloc,c),score(oloc,c));
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
newyl = [yl(1) yl(1) + 1.17*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,sprintf('p = %1.3f',p),'horizontalalignment','center','fontsize',15)
title(sprintf('Component %d score by localization',c))
set(gca,'fontsize',15)
ylim(newyl)
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})

nexttile([1 2])
iqr_spikes = prctile(orig_bins,[25 75],1);
median_spikes = nanmedian(orig_bins,1);
shaded_error_bars(times,median_spikes,iqr_spikes,[0 0 0]);
hold on
plot([0 0],ylim,'k--','LineWidth',2);
xlim([times(1) times(end)])
title('All patients')
set(gca,'fontsize',15)
xlabel('Hours')
ylabel('Spikes/elecs/min')

for i = 1:4
    nexttile
    plot(times,orig_bins(top_12(i),:),'linewidth',2)
    hold on
    plot([0 0],ylim,'k--','linewidth',2)
    if i <= 2
        top_text = 'highest';
    else
        top_text = 'lowest';
    end
    if i == 1 || i == 3
        first_text = '1st';
    else
        first_text = '2nd';
    end
    title(sprintf('Patient with %s\n%s component score',top_text,first_text))
    set(gca,'fontsize',15)
    xlabel('Hours')
    ylabel('Spikes/elecs/min')
    xlim([times(1) times(end)])
end


title(t,title_text,'fontsize',20,'fontweight','bold')
end