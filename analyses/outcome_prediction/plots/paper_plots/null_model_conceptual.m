function null_model_conceptual

which_atlas = 'brainnetome';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [locations.main_folder,'data/atlas/brainnetome/'];
atlas_out_folder = [results_folder,'analysis/atlas/'];

bct_folder= locations.bct;
data_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/paper_plots/'];

%% Load stuff
data = load([data_folder,'main_out.mat']);
data = data.out;

model_folder = [results_folder,'analysis/outcome/plots/'];
model = load([model_folder,'model_stuff_',which_atlas,'.mat']);
model_info = model.all_out.model_info;

A = niftiread([atlas_folder,'BN_Atlas_246_2mm.nii.gz']);

atlas_out = load([atlas_out_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;


p = 80;

%% Do atlas prep
atlas_nums = atlas_out.atlas_nums;
atlas_names = atlas_out.atlas_names;
broad = localize_regions(atlas_names,which_atlas);
broad_no_lat = cell(length(broad),1);
for ib = 1:length(broad)
    curr = broad{ib};
    if isempty(curr),continue; end
    curr = strrep(curr,'left ','');
    curr = strrep(curr,'right ','');
    broad_no_lat{ib} = curr;
end
broad_nums = nan(size(broad_no_lat));
broad_nums(strcmp(broad_no_lat,'mesial temporal')) = 1;
broad_nums(strcmp(broad_no_lat,'temporal neocortical')) = 2;
broad_nums(strcmp(broad_no_lat,'other cortex')) = 3;
%broad_nums(isnan(broad_nums)) = 4;

broad_A = nan(size(A));
non_zero_A = A~=0;
broad_A(non_zero_A) = broad_nums(A(non_zero_A));
broad_A(isnan(broad_A)) = 0;

%% Prep figure
fig_name = 'Fig 5';
figure
set(gcf,'position',[10 10 1000 700])
tiledlayout(2,2,'tilespacing','tight','padding','tight')


%% A: Brain showing MT, TN, Other
nexttile
s = 30;
slice = broad_A(:,:,s);
slice = imrotate(slice,90);
%slice = imresize(slice,'scale',[0.8 1]);
image(int8(slice));
cmap = colormap(lines);
cmap = [1 1 1;cmap];
colormap(gca,cmap);
axis square;
%colormap(lines)
grid off
set(gca,'visible','off')
text(size(slice,2)/2,5,'Anatomical regions','fontsize',20,...
    'horizontalalignment','center','fontweight','bold')

%% B: Conceptual figure of how to calculate density
nexttile
r = 20;
x = (-40:40)';
y = (-40:40)';
[X,Y] = meshgrid(x,y);
d = sqrt(X.^2 + Y.^2);
Z = density_function(d,r);
range = max(Z(:))-min(Z(:));
maxZ = max(Z(:));

surf(x,y,Z)
colormap(gca,parula)
hold on
plot3([0 r],[0 0],[maxZ + 0.05*range maxZ + 0.05*range],'k','linewidth',2)
plot3([0 0],[0 0],[maxZ + 0.01*range maxZ + 0.09*range],'k','linewidth',2)
plot3([r r],[0 0],[maxZ + 0.01*range maxZ + 0.09*range],'k','linewidth',2)
text(r/2,0,maxZ + 0.15*range,'Search radius','fontsize',20,...
    'horizontalalignment','center','rotation',0)
view(0,15)
xticklabels([])
yticklabels([])
zticklabels([])
grid off
set(gca,'visible','off')
text(0,0,maxZ + 0.28*range,'Spatial density model','fontsize',20,...
    'horizontalalignment','center','fontweight','bold')

%% C: Spatial density example (single patient)
locs = data.all_locs{p};
sr = calculate_default_search_radius(data.all_locs);
density = estimate_coverage_density(locs,sr);
nexttile
scatter3(locs(:,1),locs(:,2),locs(:,3),100,density,'filled',...
    'markeredgecolor','k','linewidth',2);
colormap(gca,parula)
%c = colorbar;
xticklabels([])
yticklabels([])
zticklabels([])
%ylabel(c,'Density (1/mm^2)')
set(gca,'fontsize',15)
text(mean(locs(:,1)),mean(locs(:,2)),max(locs(:,3))+45,...
    sprintf('Electrode coverage density\n(single patient example)'),'fontsize',20,...
    'horizontalalignment','center','fontweight','bold')
grid off
set(gca,'visible','off')

%% D: ROC for null model
nexttile
pp = plot(model_info(2).x,model_info(2).ym,'k','linewidth',2);
hold on
plot([0 1],[0 1],'k:','linewidth',2)
set(gca,'fontsize',20)
xlabel('False positive rate')
ylabel('True positive rate')
legend(pp,sprintf('AUC %1.2f',mean(model_info(2).all_auc)),'fontsize',20,'location','southeast')
title('Null model ROC')

annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.46 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.46 0.1 0.1],'String','D','fontsize',25,'linestyle','none')

print(gcf,[plot_folder,fig_name],'-dpng')

end


function curr_dens = density_function(d,r)

W = 1;
%if d > r
%    curr_dens = 0;
%else
%    curr_dens = 3/pi * W * (1 - (d/r).^2).^2;
%end
curr_dens = (1/r^2*3/pi * W * (1 - min(1,(d/r)).^2).^2)/r^2;
end