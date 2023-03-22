function plot_elecs_on_brain

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
icbm_folder = [locations.main_folder,'data/ICBM152_2023/'];
box_path = locations.box_folder;
elec_path = [box_path,'CNT Implant Reconstructions/'];
pial_folder = [locations.main_folder,'data/example_surf/RID405/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Read pial files
lobj = SurfStatReadSurf([pial_folder,'lh.pial']);
lvertices = lobj.coord';
lfaces = lobj.tri;

robj = SurfStatReadSurf([pial_folder,'rh.pial']);
rvertices = robj.coord';
rfaces = robj.tri;

% try to combine left and right
vertices = [lvertices;rvertices];
faces = [lfaces;rfaces];

% Load the elecs
%{
T = readtable([pt_folder,'electrodenames_coordinates_native.csv']);
locs = [T.Var4 T.Var2 T.Var3];
names = T.Var1;
%}
%T = readtable([pial_folder,'sub-RID0405_ses-clinical01_space-T00mri_atlas-atropos_radius-2_desc-vox_coordinates.csv']);
T = readtable([pial_folder,'sub-RID0405_ses-clinical01_space-T00mri_atlas-DKTantspynet_radius-2_desc-vox_coordinates.csv']);
names = T.name;
offset = [115,-130,-80];
locs = [-T.x T.y T.z];
locs = locs + repmat(offset,size(locs,1),1);
allowable_labels = get_allowable_elecs('HUP100');
mt_symm = find_mt_symmetric_coverage(names,allowable_labels);
mt = ismember(names,mt_symm);


figure
rh = trisurf(lfaces,lvertices(:,1),lvertices(:,2),lvertices(:,3));
hold on
lh = trisurf(rfaces,rvertices(:,1),rvertices(:,2),rvertices(:,3));
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
            toff = [15 0 0];
        else
            toff = [-15 0 0];
        end
        curr_loc = curr_loc+toff;
        
        if ismember(curr,mt_symm)
            text(curr_loc(1),curr_loc(2),curr_loc(3),curr(1:2),'HorizontalAlignment','center','fontsize',25,'color','r')
        else
            text(curr_loc(1),curr_loc(2),curr_loc(3),curr(1:2),'HorizontalAlignment','center','fontsize',25,'color','k')
        end
        
    end
end
%text(locs(:,1),locs(:,2),locs(:,3),names,'horizontalalignment','center')
view(2.4,-3.6)
axis off
% Show them together
%{
figure
rh = trisurf(lfaces,lvertices(:,1),lvertices(:,2),lvertices(:,3));
hold on
lh = trisurf(rfaces,rvertices(:,1),rvertices(:,2),rvertices(:,3));
rh.LineStyle = 'none';
rh.FaceAlpha = 0.2;
rh.FaceColor = [1 0.7 0.7];    

lh.LineStyle = 'none';
lh.FaceAlpha = 0.2;
lh.FaceColor = [0.7 1 0.7];  
hold on
scatter3(locs(right_elecs,1),locs(right_elecs,2),locs(right_elecs,3),'MarkerFaceColor','w','MarkerEdgeColor','r');
scatter3(locs(~right_elecs,1),locs(~right_elecs,2),locs(~right_elecs,3),'MarkerFaceColor','w','MarkerEdgeColor','g');

text(locs(:,1),locs(:,2),locs(:,3),names,'horizontalalignment','center')
%}

end