function epilepsia_supplemental_fig1

fsize = 20;
myColours = [0 33 87;...
122 80 113;...    
227 124 29;...
    86 152 163]/255;



locations = fc_toolbox_locs;
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/epilepsia/'];

X = seeg_ad_out.X;
Y = seeg_ad_out.Y;
altX = seeg_ad_out.altX;
altY = seeg_ad_out.altY;
AUC = seeg_ad_out.AUC;
altAUC = seeg_ad_out.altAUC;
scores = seeg_ad_out.scores;
all_out = seeg_ad_out.all_out;
all_ss = seeg_out.all_ss;


% Get identities of each sleep stage
ss = all_out(:,2);

nss = length(all_ss);
ss_idx = nan(nss,length(ss));
for is = 1:nss
    ss_idx(is,:) = strcmp(ss,all_ss{is});
end


figure
set(gcf,'position',[289 517 1350 400])
tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

% Histogram of distribution of states
nexttile
b = bar(sum(ss_idx,2));
b.FaceColor = 'flat';
b.CData = repmat([0 118 192]/255,5,1);
b.CData(3,:) = [163 2 52]/255;
xticklabels(all_ss)
ylabel('Number of segments')
xlabel('SleepSEEG Sleep stage')
title({'Number of segments','in each sleep stage'})
set(gca,'fontsize',fsize)


% ROC of sleep vs wake by AD ratio
nexttile
plot(X,Y,'linewidth',2,'color',[0 118 192]/255)
hold on
plot(altX,altY,'linewidth',2,'Color',[103 119 25]/255)
xlabel('False positive rate')
ylabel('True positive rate')
set(gca,'fontsize',fsize)
title({'ROC classifying','SleepSEEG sleep by normalized ADR'})
legend({sprintf('All sleep vs wake: AUC %1.2f',AUC),...
    sprintf('N3 vs wake: AUC %1.2f',altAUC)},'location','southeast',...
    'fontsize',fsize)

% AD ratio of different stage
nexttile
b=boxplot(scores(~strcmp(ss,'N1')),ss(~strcmp(ss,'N1')),'grouporder',...
    {'R','W','N2','N3'},'Colors',myColours);
set(b,'linewidth',2)
set(gca,'fontsize',fsize)
ylabel('Normalized ADR')
xlabel('SleepSEEG classification')
h = findobj(gcf,'tag','Outliers');
for i = 1:numel(h)
    xpos = h(i).XData(1);
    if isnan(xpos), continue; end
    h(i).MarkerEdgeColor = myColours(xpos,:);
end

title({'Normalized ADR','by SleepSEEG classification'})
fontname(gcf,"calibri");
print(gcf,[out_folder,'FigS1'],'-depsc')
print(gcf,[out_folder,'FigS1'],'-dpng')

close all
end