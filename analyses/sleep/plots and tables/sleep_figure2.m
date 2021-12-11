function sleep_figure2

%% Parameters
plot_type = 'scatter';
nblocks = 6;
myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];



locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

%% Load out file and get roc stuff
out = load([out_folder,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

figure
set(gcf,'position',[100 100 900 700])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

%% Sleep surge histogram
all_pts_spikes_bins = sleep_hist_out.all_pts_spikes_bins;
all_pts_sleep_bins = sleep_hist_out.all_pts_sleep_bins;
sig_bins = sleep_hist_out.sig_bins;

ax1 = nexttile;
median_spikes = nanmedian(all_pts_spikes_bins,1);
mean_sleep = nanmean(all_pts_sleep_bins,1)*100;
iqr_spikes = prctile(all_pts_spikes_bins,[25,75],1);
times = linspace(-12,12,length(median_spikes));
ylabels = ["Spikes/elecs/min","% Asleep"];
s = stackedplot(times,[median_spikes',mean_sleep'],'linewidth',2,...
    "DisplayLabels",ylabels);

for k = 1:length(s.LineProperties)
    if k <= size(myColours,1)
        s.LineProperties(k).Color = myColours(k,:);
    end
end
ax = findobj(s.NodeChildren, 'Type','Axes');
arrayfun(@(h)xline(h,0,'--k','LineWidth',2),ax)
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
xlabel('Hours relative to sleep onset')
xlim([-12 12])
set(gca,'fontsize',15)
title('Spike rate and proportion asleep around sleep onset')

%% sleep pca stuff
sleep_bins = out.sleep_hist_out.all_pts_spikes_bins;
% Subtract mean
sleep_bins = (sleep_bins - nanmean(sleep_bins,2))./nanstd(sleep_bins,[],2);
[coeff,score,latent] = pca(sleep_bins,'Rows','pairwise');
locs = out.circ_out.all_locs;
% Top scorers
[~,top] = max(score,[],1);
% Bottom scorers
[~,bottom] = min(score,[],1);

if 1
    f = figure;
    turn_nans_gray(sleep_bins)
    xlabel('Time bins')
    ylabel('Patients')
    title('Normalized spike rates')
    colorbar
    set(gca,'fontsize',15)
    print(f,[out_folder,'pca_fig'],'-dpng')
    close(gcf)
end

%% Variances
%
ax2 = nexttile;
stem(latent,'k','linewidth',2);
%}


%% First two components
ax3 = nexttile;
ylabels = ["Component 1","Component 2"];
sh = stackedplot(times,[coeff(:,1),coeff(:,2)],'linewidth',2,...
    "DisplayLabels",ylabels);
for k = 1:length(sh.LineProperties)
    if k <= size(myColours,1)
        sh.LineProperties(k).Color = myColours(k+2,:);
    end
end
axh = findobj(sh.NodeChildren, 'Type','Axes');
arrayfun(@(h)xline(h,0,'--k','LineWidth',2),axh)
%pause(0.5)
set([axh.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
xlabel('Hours relative to sleep onset')
xlim([-12 12])
set(gca,'fontsize',15)
title('Top 2 principal component coefficients')
%}

%% F: Localization
ax4 = nexttile;
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
text(1.5,ytext,get_p_text(p),'horizontalalignment','center','fontsize',15)
title(sprintf('Component %d score by localization',c))
set(gca,'fontsize',15)
ylim(newyl)
xticks([1 2])
ylabel('Score')
xticklabels({'Temporal','Extra-temporal'})

%% Fix for wacky error - add labels and such later
title(ax2,'Principal component variances')
set(ax2,'fontsize',15)
xlabel(ax2,'Component')
ylabel(ax2,'Variance')

%% Correlate other things with PCA score
sex = circ_out.sex;
age_onset = circ_out.age_onset;
age_implant = circ_out.age_implant;
duration = circ_out.duration;

% Sex
male = strcmp(sex,'Male');
female = strcmp(sex,'Female');


print([out_folder,'Fig2'],'-dpng')

end