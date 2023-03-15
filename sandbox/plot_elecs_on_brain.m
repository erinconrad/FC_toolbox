function plot_elecs_on_brain

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
icbm_folder = [locations.main_folder,'data/ICBM152_2023/'];
box_path = locations.box_folder;
elec_path = [box_path,'CNT Implant Reconstructions/'];
pt_folder = [elec_path,'RID405_HUP190/'];
mprage_file = [pt_folder,'T00_RID405_mprage_brainBrainExtractionBrain.nii.gz'];
be_file = [pt_folder,'T00_RID405_mprage_brainBrainExtractionBrain/T00_RID405_mprage_brainBrainExtractionBrain_wholebrainseg.nii.gz'];
pial_folder = [locations.main_folder,'data/example_surf/RID405/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Read pial files
lobj = readObj([pial_folder,'lh.obj']);
lvertices = lobj.v;
lfaces = lobj.f.v;

robj = readObj([pial_folder,'rh.obj']);
rvertices = robj.v;
rfaces = robj.f.v;

% try to combine left and right
vertices = [lvertices;rvertices];
faces = [lfaces;rfaces];

% Load the elecs
%{
T = readtable([pt_folder,'electrodenames_coordinates_native.csv']);
locs = [T.Var4 T.Var2 T.Var3];
names = T.Var1;
%}
T = readtable([pial_folder,'sub-RID0405_ses-clinical01_space-T00mri_atlas-atropos_radius-2_desc-vox_coordinates.csv']);
names = T.name;
locs = [T.x T.y T.z];
right_elecs = cellfun(@(x) strcmp(x(1),'R'),names);


% Show them together
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

end