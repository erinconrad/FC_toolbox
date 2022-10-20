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

figure
set(gcf,'position',[289 517 1001 280])
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile
plot(X,Y,'linewidth',2)
hold on
plot(altX,altY,'linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
set(gca,'fontsize',15)
title('ROC classifying SleepSEEG sleep by normalized ADR')
legend({sprintf('All sleep vs wake AUC %1.2f',AUC),sprintf('N3 vs wake AUC %1.2f',altAUC)},'location','southeast')


nexttile
boxplot(scores,all_out(:,2))
set(gca,'fontsize',15)
ylabel('Normalized alpha-delta ratio')
xlabel('SleepSEEG classification')
title('Normalized ADR by SleepSEEG classification')



end