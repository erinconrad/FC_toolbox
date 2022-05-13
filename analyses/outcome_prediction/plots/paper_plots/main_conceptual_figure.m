function main_conceptual_figure

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
bct_folder= locations.bct;
atlas_folder = [results_folder,'analysis/atlas/'];
data_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/paper_plots/'];
model_folder = [results_folder,'analysis/outcome/plots/'];

if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load stuff
data = load([data_folder,'main_out.mat']);
data = data.out;

%% initialize figure
figure
set(gcf,'position',[10 10 1000 700])
t1 = tiledlayout(4,2,'tilespacing','compact','padding','tight');


%% A: incomplete electrode coverage
offset = [0 5 35];
name = 'HUP106';
brainFolder = [locations.main_folder,'data/pretty_brain/'];
giftiFolder = [brainFolder,name,'/'];
names = dir([giftiFolder,'*pial']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);
ip = strcmp(data.all_names,name);
%T = readtable([giftiFolder,'electrode_coordinates_native.csv']);
%locs = [T.Var1 T.Var2 T.Var3];
locs = data.all_locs{ip};

locs(any(locs>1e5,2),:) = nan;
locs = locs + repmat(offset,size(locs,1),1);

ax3 = nexttile([2 1]);
plot_gifti_erin_2(ax3,g)

hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'filled','markerfacecolor','k',...
    'markeredgecolor','k')
view([-99 7])


%% B: normal anatomical variation in connectivity
which_atlas = 'aal_bernabei';
rate = data.all_spikes;
soz = data.all_soz_bin;
npts = length(soz);
labels = data.all_labels;
fc = data.all_fc;
coh = data.all_coh;
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;
atlas_nums = atlas_out.atlas_nums;
names = atlas_out.atlas_names;
[locs,lats] = lateralize_regions(names,which_atlas);
atlas = rebuild_atlas(fc,rate,atlas_elec_labels,...
    atlas_elec_regions,atlas_nums,labels,soz,coh);
lr_order = reorder_lr(locs,lats);
atlas = atlas(lr_order,lr_order,:);

nexttile([2 1])
mean_atlas = nanmean(atlas,3);
nan_regions = sum(~isnan(mean_atlas),1) == 0;
mean_atlas(nan_regions,:) = [];
mean_atlas(:,nan_regions) = [];
turn_nans_gray(mean_atlas)
xticklabels([])
yticklabels([])
c = colorbar('location','eastoutside');
ylabel(c,'Correlation (r)')
%xlabel('Anatomical region')
ylabel('Anatomical region')
title('Anatomical variability in connectivity')
set(gca,'fontsize',20)

%% C: Density
p = 102;
locs = data.all_locs{p};
sr = calculate_default_search_radius(data.all_locs);
density = estimate_coverage_density(locs,sr);
ax4 = nexttile([2 1]);
scatter3(locs(:,1),locs(:,2),locs(:,3),100,density,'filled',...
    'markeredgecolor','k','linewidth',2);
colormap(gca,parula)
c = colorbar('location','southoutside');
xticklabels([])
yticklabels([])
zticklabels([])
ylabel(c,'Density (1/mm^2)')
set(gca,'fontsize',20)

grid off
set(gca,'visible','off')
view(-57.5,-23)

%% D: SPikes
spout = get_spikes_and_corr;
times = spout.times;
sp_chs = spout.sp_chs;
values = spout.values;
mean_times = spout.mean_times;
all_corrs = spout.all_corrs;
window = spout.window;


nexttile
offset = 0;
for isp = 1:length(sp_chs)
    plot(times,values(:,sp_chs(isp))-offset,'k','linewidth',2)
    hold on
    if isp < length(sp_chs)
        offset = offset - (min(values(:,sp_chs(isp))) - max(values(:,sp_chs(isp+1))));
    end
end
title('Spikes alter connectivity')
xticklabels([])
yticklabels([])
ylabel('Amplitude')
set(gca,'fontsize',20)

nexttile
plot(mean_times,all_corrs,'linewidth',3)
xlim([0 window])
yticklabels([])
xlabel('Time (s)')
ylabel('Correlation (r)')
set(gca,'fontsize',20)

%% text first plot
text(ax3,mean(locs(:,1)),mean(locs(:,2)),max(locs(:,3))+92,...
    sprintf('Incomplete electrode sampling'),'fontsize',23,...
    'horizontalalignment','center','fontweight','bold')
text(ax4,mean(locs(:,1)),mean(locs(:,2))-15,max(locs(:,3))+49,...
    sprintf('Variability in electrode density'),'fontsize',23,...
    'horizontalalignment','center','fontweight','bold')

%% Annotations
annotation('textbox',[0 0.90 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.47 0.90 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.46 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.47 0.46 0.1 0.1],'String','D','fontsize',25,'linestyle','none')

print(gcf,[plot_folder,'Fig 1'],'-dpng')

end

