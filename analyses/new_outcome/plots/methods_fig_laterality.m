function methods_fig_laterality

%{
Plan:
A: Show a brain with electrode coordinates, highlighting the mesial
temporal targets
B: Show snippets of EEG, with arrows going to a list of features being
extracted, showing sleep stages, references, frequency bands, average + SD
C: Show calculation of asymmetry index
D: Show example of AI according to epilepsy laterality
%}

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/Methods/'];
pial_folder = [locations.main_folder,'data/example_surf/RID652/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
freesurfer_path = '/Applications/freesurfer/7.3.2/matlab/';

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

addpath(genpath(freesurfer_path))

if 1
%% Read pial files

lobj = SurfStatReadSurf([pial_folder,'lh.pial']);
lvertices = lobj.coord';
lfaces = (lobj.tri);

robj = SurfStatReadSurf([pial_folder,'rh.pial']);
rvertices = robj.coord';
rfaces = robj.tri;

%% Get electrodes
%T1 = readtable([pial_folder,'sub-RID0296_ses-clinical01_space-T00mri_atlas-DKTantspynet_radius-2_desc-vox_coordinates.csv']);
%names = T1.name;
T1 = readtable([pial_folder,'sub-RID0652_electrode_names.txt'],'ReadVariableNames',false);
names = T1.Var1;

%T = readtable([pial_folder,'sub-RID0652_ses-clinical01_space-T00mri_desc-mm_electrodes.txt']);
T = readtable([pial_folder,'sub-RID0652_ses-research3T_space-T00mri_desc-mm_electrodes.txt']);

locs = [T.Var1 T.Var2 T.Var3];
allowable_labels = get_allowable_elecs('HUP100');
mt_symm = find_mt_symmetric_coverage(names,allowable_labels);
mt = ismember(names,mt_symm);

%% Load the T1 to get information to convert pial surface to electrode space
mri = MRIread([pial_folder,'T1.mgz']);
vox_2_ras = mri.vox2ras;
tkras = mri.tkrvox2ras;
transform = @(x) ((vox_2_ras * inv(tkras) * ([x repmat(ones,size(x,1),1)])')');
first_three_columns = @(x) x(:,1:3);
all_trans = @(x) first_three_columns(transform(x));

rvt = all_trans(rvertices);
lvt = all_trans(lvertices);



figure
rh = trisurf(lfaces,lvt(:,1),lvt(:,2),lvt(:,3));
hold on
lh = trisurf(rfaces,rvt(:,1),rvt(:,2),rvt(:,3));
hold on
rh.LineStyle = 'none';
rh.FaceAlpha = 0.1;
rh.FaceColor = [0.7 0.6 0.6];
lh.LineStyle = 'none';
lh.FaceAlpha = 0.1;
lh.FaceColor = [0.7 0.6 0.6];
scatter3(locs(~mt,1),locs(~mt,2),locs(~mt,3),'markerfacecolor','k','markeredgecolor','k')
scatter3(locs(mt,1),locs(mt,2),locs(mt,3),'markerfacecolor','r','markeredgecolor','r')


% Name LA, LB, etc.
end_elecs = {'LA12','LB12','LC12','RA12','RB12','RC12'};
for i = 1:length(end_elecs)
    curr = end_elecs{i};
    match = strcmp(names,curr);
    if sum(match) ~=0
        curr_loc = locs(match,:);
        if strcmp(curr(1),'L')
            toff = [-10 0 0];
        else
            toff = [10 0 0];
        end
        curr_loc = curr_loc+toff;
        
        if ismember(curr,mt_symm)
            text(curr_loc(1),curr_loc(2),curr_loc(3),curr(1:2),'HorizontalAlignment','center','fontsize',25,'color','r')
        else
            text(curr_loc(1),curr_loc(2),curr_loc(3),curr(1:2),'HorizontalAlignment','center','fontsize',25,'color','k')
        end
        
    end
end

view(-180,-90)%view(-176,4) %view(-182,-5)
axis off
print(gcf,[plot_folder,'elecs_example'],'-dpng')
print(gcf,[plot_folder,'elecs_example'],'-depsc')

close(gcf)
end

%% B snippet of eeg with features

%% Asymmetry index

%% Now do lr_mt to get AI features
[T,features] =  lr_mt; 
response = 'soz_lats';
figure
%% Show spikes
feature = 'spikes car sleep';
h = boxplot_with_points(T.(feature),T.(response),0,{'left','right','bilateral'});
spikes = T.(feature);
left_soz_spikes = spikes(strcmp(T.(response),'left'));
right_soz_spikes = spikes(strcmp(T.(response),'right'));
ylabel('Spike rate asymmetry index')
%title('Spike rate asymmetry index by SOZ laterality')
set(gca,'fontsize',20)
print(gcf,[plot_folder,'AI_example'],'-dpng')
print(gcf,[plot_folder,'AI_example'],'-depsc')

end