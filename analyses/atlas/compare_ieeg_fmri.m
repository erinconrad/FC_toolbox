function compare_ieeg_fmri

which_atlas = 'brainnetome';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];
bct_folder= locations.bct;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load ieeg atlas
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;

atlas = out.atlas;
names = out.atlas_names;

if size(names,1) < size(names,2)
    names = names';
end
ieeg = out.atlas;
all_ieeg = ieeg;
ieeg = nanmean(ieeg,3);

%% Measure sparsity
all_ieegw = wrap_or_unwrap_adjacency_fc_toolbox(all_ieeg);
pts_per_edge = measure_sparsity(all_ieegw);
if 0
    histogram(pts_per_edge)
end

%% Load fmri atlases
fmri = load([atlas_folder,which_atlas,'_sc_fc.mat']);
fmri = fmri.out;
fmri_fc = (fmri.FC);
fmri_sc = fmri.SC;

%% Find and remove nan rows and things not usually targeted with implantation
basal_ganglia = contains(names,'BG_');
thalamus = contains(names,'Tha_');
nan_rows = (sum(isnan(ieeg),1) == size(ieeg,2))';
remove = basal_ganglia | thalamus | nan_rows;

fmri_fc = fmri_fc(~remove,~remove);
fmri_sc = fmri_sc(~remove,~remove);
ieeg = ieeg(~remove,~remove);

%% Do corr
ieegw = wrap_or_unwrap_adjacency_fc_toolbox(ieeg);
fmri_fcw = wrap_or_unwrap_adjacency_fc_toolbox(fmri_fc);
fmri_scw = wrap_or_unwrap_adjacency_fc_toolbox(fmri_sc);
[r_fc,p_fc] = corr(ieegw,fmri_fcw,'rows','pairwise');
[r_sc,p_sc] = corr(ieegw,fmri_scw,'rows','pairwise');

%% Show side by side and corr
figure
set(gcf,'position',[10 10 900 900])
tiledlayout(2,3,'tilespacing','tight','padding','tight')

nexttile
imagesc(fmri_fc)
title('fMRI FC')

nexttile
imagesc(ieeg)
title('IEEG')

nexttile
plot(ieegw,fmri_fcw,'o')
xl = xlim;
yl = ylim;
text(xl(2),yl(1),sprintf('r = %1.2f, %s',r_fc,get_p_text(p_fc)),...
    'horizontalalignment','right','verticalalignment','bottom');
title('IEEG-fMRI correlation')

nexttile
imagesc(fmri_sc)
title('dMRI')

nexttile
imagesc(ieeg)
title('IEEG')

nexttile
plot(ieegw,fmri_scw,'o')
xl = xlim;
yl = ylim;
text(xl(2),yl(1),sprintf('r = %1.2f, %s',r_sc,get_p_text(p_sc)),...
    'horizontalalignment','right','verticalalignment','bottom');
title('IEEG-dMRI correlation')

end