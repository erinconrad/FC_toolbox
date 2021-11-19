function any_pca(spike_bins,times,locs,names)

%ex_pts = [21,55,23,13,50,77];

npts = size(spike_bins,1);
nbins = size(spike_bins,2);

% Subtract mean
spike_bins = (spike_bins - nanmean(spike_bins,2))./nanstd(spike_bins,[],2);

[coeff,score,latent] = pca(spike_bins,'Rows','complete');

% Top scorers
[~,top] = max(score,[],1);

% Bottom scorers
[~,bottom] = min(score,[],1);

% rm top scorers
%{
allpts = 1:npts;
rm_top = allpts;
rm_top(ismember(rm_top,t1(1:3))) = [];
[~,t2] = max(score(rm_top,:),[],1);
t2 = rm_top(t2);
%}


top_12 = [top(1:3),bottom(1:3)];

figure
set(gcf,'position',[10 10 1100 400])
t = tiledlayout(2,6,'tilespacing','tight','padding','tight');

for i = 1:6
    nexttile
    plot(times,spike_bins(top_12(i),:),'linewidth',2)
    hold on
    plot([0 0],ylim,'k--','linewidth',2)
    title(names(top_12(i)))
end

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
p3 = plot(times,coeff(:,3),'linewidth',2);
title('Top principal components')
legend([p1 p2 p3],{'Component 1','Component 2','Component 3'},'fontsize',15,'location','northwest')
set(gca,'fontsize',15)
xlabel('Hours')

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
ytext = yl(1) + 1.09*(yl(2)-yl(1));
newyl = [yl(1) yl(1) + 1.15*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,sprintf('p = %1.3f',p),'horizontalalignment','center','fontsize',15)
title(sprintf('Component %d score by localization',c))
set(gca,'fontsize',15)
ylim(newyl)
xticks([1 2])
xticklabels({'Temporal','Other'})


title(t,'Peri-ictal spike patterns','fontsize',20,'fontweight','bold')
end