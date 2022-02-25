function laterality_fig(in)

plot_type = 'scatter';

unpack_any_struct(in);

figure
set(gcf,'position',[10 10 1100 1100])
tiledlayout(3,3,'tilespacing','tight','padding','tight')

%{
%% L-R ordered atlas
nexttile
thing = nanmean(atlas,3);
all_nan = sum(isnan(thing),2)==size(thing,1);
thing(all_nan,:) = []; thing(:,all_nan) = [];
left = zeros(size(atlas,1),1);
left(1:size(atlas,1)/2) = 1;
left(all_nan) = [];
cutoff = find(left);
cutoff = cutoff(end);
%{
temp_names = names; temp_names(all_nan) = [];
table(temp_names,(left))
%}
pretty_matrix(thing,{'Left','Right'},cutoff,'r^2')
title('Average functional connectivity')
set(gca,'fontsize',15)
%}

%% Left vs right intrinisic connectivity
nexttile
stats = plot_paired_data((bl_intrinsic)',{'left','right','right'},'Intrinsic connectivity','paired',plot_type);
title('Left vs right intrinsic connectivity')
set(gca,'fontsize',15)


%% SOZ-non SOZ ordered atlas
nexttile
thing = nanmean(soz_atlas,3);
all_nan = sum(isnan(thing),2)==size(thing,1);
thing(all_nan,:) = []; thing(:,all_nan) = [];
soz = zeros(size(soz_atlas,1),1);
soz(1:size(soz_atlas,1)/2) = 1;
soz(all_nan) = [];
cutoff = find(soz);
cutoff = cutoff(end);
pretty_matrix(thing,{'SOZ','Non-SOZ'},cutoff,'r^2')
title('Average functional connectivity')
set(gca,'fontsize',15)

%%  SOZ vs on SOZ instrinsic connectivity
nexttile
stats = plot_paired_data((soz_non_intrinsic)',{'SOZ','non-SOZ','non-SOZ'},'Connectivity','paired',plot_type);
title('SOZ vs non-SOZ intrinsic connectivity')
set(gca,'fontsize',15)

%% SOZ-non SOZ coverage map
nexttile
thing = cov_map;
all_nan = sum(isnan(thing),2)==size(thing,1);
thing(all_nan,:) = []; thing(:,all_nan) = [];
soz = zeros(size(soz_atlas,1),1);
soz(1:size(soz_atlas,1)/2) = 1;
soz(all_nan) = [];
cutoff = find(soz);
cutoff = cutoff(end);
pretty_matrix(thing,{'SOZ','Non-SOZ'},cutoff,'# Patients with coverage')
title('Patient electrode coverage')
set(gca,'fontsize',15)



%% SOZ-non SOZ ordered atlas symmetric
nexttile
thing = nanmean(symm_soz_atlas,3);
all_nan = sum(isnan(thing),2)==size(thing,1);
thing(all_nan,:) = []; thing(:,all_nan) = [];
soz = zeros(size(soz_atlas,1),1);
soz(1:size(soz_atlas,1)/2) = 1;
soz(all_nan) = [];
cutoff = find(soz);
cutoff = cutoff(end);
pretty_matrix(thing,{'SOZ','Non-SOZ'},cutoff,'r^2')
title({'Average functional connectivity','(symmetric electrode coverage)'})
set(gca,'fontsize',15)

%{
%% SOZ-non SOZ coverage map symmetric
nexttile
thing = symm_cov;
all_nan = sum(isnan(thing),2)==size(thing,1);
thing(all_nan,:) = []; thing(:,all_nan) = [];
soz = zeros(size(soz_atlas,1),1);
soz(1:size(soz_atlas,1)/2) = 1;
soz(all_nan) = [];
cutoff = find(soz);
cutoff = cutoff(end);
pretty_matrix(thing,{'SOZ','Non-SOZ'},cutoff,'# Patients with symmetric coverage')
title('Symmetric electrode coverage')
set(gca,'fontsize',15)
%}

%%  SOZ vs on SOZ instrinsic connectivity
nexttile
stats = plot_paired_data((symm_soz_not)',{'SOZ','non-SOZ','non-SOZ'},'Connectivity','paired',plot_type);
title({'SOZ vs non-SOZ intrinsic connectivity', '(symmetric coverage)'})
set(gca,'fontsize',15)

%% FC confusion matrix
nexttile
plot_confusion_matrix(fc_conf,{'Right','Left'},'SOZ predicted by low connecitivity','True SOZ')
title('Connectivity accuracy at predicting SOZ')
set(gca,'fontsize',15)

%% FC confusion matrix
nexttile
plot_confusion_matrix(spikes_conf,{'Right','Left'},'SOZ predicted by frequent spikes','True SOZ')
title('Spike accuracy at predicting SOZ')
set(gca,'fontsize',15)

%% Something predicting bilaterality
nexttile
plot(1+randn(sum(unilat),1)*0.05,Q(unilat),'o','linewidth',2)
hold on
plot(2+randn(sum(bilat),1)*0.05,Q(bilat),'o','linewidth',2)
xticks([1 2])
xlim([0.5 2.5])
xticklabels({'Unilateral','Bilateral'})
ylabel('Modularity')
title('Network modularity by laterality')
set(gca,'fontsize',15)
[p,~,stats] = ranksum(Q(unilat),Q(bilat));

yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)

print(gcf,[plot_folder,'fig1'],'-dpng')


end