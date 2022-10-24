function epilepsia_supplemental_fig1

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
set(gcf,'position',[289 517 1350 280])
tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

% Histogram of distribution of states
nexttile
b = bar(sum(ss_idx,2));
b.FaceColor = 'flat';
b.CData(3,:) = [1 0 0];
xticklabels(all_ss)
ylabel('Number of segments')
xlabel('SleepSEEG Sleep stage')
title('Number of segments in each sleep stage')
set(gca,'fontsize',15)


% ROC of sleep vs wake by AD ratio
nexttile
plot(X,Y,'linewidth',2)
hold on
plot(altX,altY,'linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
set(gca,'fontsize',15)
title('ROC classifying SleepSEEG sleep by normalized ADR')
legend({sprintf('All sleep vs wake: AUC %1.2f',AUC),sprintf('N3 vs wake: AUC %1.2f',altAUC)},'location','southeast')

% AD ratio of different stage
nexttile
b=boxplot(scores(~strcmp(ss,'N1')),ss(~strcmp(ss,'N1')),'grouporder',...
    {'R','W','N2','N3'});
set(b,'linewidth',2)
set(gca,'fontsize',15)
ylabel('Normalized alpha-delta ratio')
xlabel('SleepSEEG classification')
title('Normalized ADR by SleepSEEG classification')

print(gcf,[out_folder,'FigS1'],'-depsc')
close all
end