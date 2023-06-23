function reproduce_problem

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
subplot_path = [plot_folder,'ai_subplots/'];
if ~exist(subplot_path,'dir')
    mkdir(subplot_path)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

figure
set(gcf,'position',[25 235 900 900])
t = tiledlayout(4,2,"TileSpacing",'tight','padding','tight');

%% fmri locs
file_path = '/Users/erinconrad/Desktop/research/FC_toolbox/Alfredo_code/fmri_analysis_AL_3_28_23/';
csv_path = [file_path,'out_csvs/'];

%% Load files
T = readtable([file_path,'df.csv']);
bT = readtable([file_path,'BNA_subregions.xlsx']);
mt_mask = niftiread([file_path,'mesial_temporal_roi.nii.gz']);
mni_brain = niftiread([file_path,'tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz']);

%% Brains
special_color = [0.9 0.1 0.1];
% make it a double
mni_brain = double(mni_brain);

% resize mni brain to match the mt mask
mni_brain = imresize3(mni_brain,size(mt_mask));

% set 0 to be the brightest
mni_brain(mni_brain <= 100) = max(mni_brain,[],'all');

% Set the mt_mask to a special value
mni_brain(mt_mask==1) = nan;

% Plots
t1 = nexttile(t,5);
colormap(t1,'gray')
turn_nans_gray(imrotate(mni_brain(:,:,60),90),special_color,t1)
axis off
title({'Regions included in fMRI','connectivity calculations'})
set(gca,'fontsize',15)

% Plot brain again
t2 = nexttile(t,7);
colormap(t2,'gray')
turn_nans_gray(imrotate(squeeze(mni_brain(:,100,:)),90),special_color,t2)
axis off



end