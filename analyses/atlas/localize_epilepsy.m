function localize_epilepsy

%% Parameters
which_atlas = 'aal_bernabei';%'brainnetome';% %'aal';'aal_bernabei';
plot_type = 'scatter';
broad_regions = {'left mesial temporal','right mesial temporal',...
    'left temporal neocortical','right temporal neocortical',...
    'left other cortex','right other cortex'};
broad_non_lat = {'mesial temporal','temporal neocortical','other cortex'};
broad_lat = {'left','right','left','right','left','right'};
nb = length(broad_regions);

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];
bct_folder= locations.bct;
out_folder = [results_folder,'analysis/atlas/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load atlas
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;

atlas = out.atlas;
names = out.atlas_names;
pt_names = out.pt_names;
atlas_nums = out.atlas_nums;
nregions = length(names);
assert(nregions==size(atlas,1))
npts = size(atlas,3);
nelecs = out.n_elecs_all;

if strcmp(which_atlas,'aal_bernabei'), names = names'; end

%% Get soz atlas identifiers
%sozs = out.sozs;

%% Get soz loc-lat combos
soz_out = load('out.mat');
soz_out = soz_out.out.circ_out;
soz_locs = soz_out.all_locs;
soz_lats = soz_out.all_lats;
soz_lat_loc = cellfun(@(x,y) [x,' ',y],soz_lats,soz_locs,'UniformOutput',false);
if 0
    table(soz_locs,soz_lats)
end


if 0
%% Get electrodes to ignore
%{
ignore = regions_to_ignore(which_atlas,names);

atlas(ignore,:,:) = nan;
atlas(:,ignore,:) = nan;
%}
end

%% Localize regions into broad categories
broad = localize_regions(names,which_atlas);
if 0
    table(names,broad)
end

%% Ignore other
%{
ignore = cellfun(@isempty,broad);
atlas(ignore,:,:) = nan;
atlas(:,ignore,:) = nan;
%}

%% Get the average connectivity of each of the broad regions to each other 
% also getting SOZ and number of electrodes
broad_connectivity = nan(nb,nb,npts);
soz_broad = zeros(nb,npts);
nelecs_broad = zeros(nb,npts);
% Loop over broad regions (6 total)
for ib = 1:nb

    % Loop over broad regions again
    for jb = 1:nb
        
        % Loop over patients
        for ip = 1:npts

            % Get the regions in the ib and jb broad regions
            i_regions = strcmp(broad,broad_regions{ib});
            j_regions = strcmp(broad,broad_regions{jb});

            % get intrinsic connectivity of that region
            conn = nanmean(atlas(i_regions,j_regions,ip),'all');
            broad_connectivity(ib,jb,ip) = conn;

            % Get soz
            curr_soz = soz_lat_loc{ip};
            soz_broad(ib,ip) = strcmp(curr_soz,broad_regions{ib});
            
            % get nelecs
            curr_nelecs = nelecs(i_regions,ip);
            nelecs_broad(ib,ip) = sum(curr_nelecs);

        end
        
    
    end
end
assert(sum(sum(soz_broad,1)<=1) == length(soz_broad))
soz_broad = logical(soz_broad);

%% Print some stats on how many SOZs in each category
for ib = 1:nb
    curr_text = broad_regions{ib};
    num_soz = sum(soz_broad(ib,:));
    fprintf('\n%d of %d (%1.1f%%) patients had SOZ in %s.\n',num_soz,npts,...
        num_soz/npts*100,curr_text);
end
fprintf('\n%d of %d (%1.1f%%) patients had SOZ that didn''t fit in above categories.\n',...
    sum(sum(soz_broad,1)==0),npts,sum(sum(soz_broad,1)==0)/npts*100);

%% Get normalized regional connectivity
% For each element of the connectivity matrix and each patient, I ask how
% strong that edge is relative to the patient average
zbroad = (broad_connectivity - nanmean(broad_connectivity,3))./nanstd(broad_connectivity,[],3);



%% Reorder connectivity matrix by SOZ - non SOZ ipsi 1 - non SOZ ipsi 2 - contralateral for each
znew = reorder_locs_soz(zbroad,soz_broad);
if 1
figure
turn_nans_gray(nanmean(znew,3))
xticks(1:6)
xticklabels({'SOZ','ipsi1','ipsi2','contra1','contra2','contra3'})
yticks(1:6)
yticklabels({'SOZ','ipsi1','ipsi2','contra1','contra2','contra3'})

end


%% Lateralize and localize epilepsy
soz_lat_not = [nanmean(squeeze([znew(1,1,:),znew(2,2,:),znew(3,3,:)]),1)',...
    nanmean(squeeze([znew(4,4,:),znew(5,5,:),znew(6,6,:)]),1)'];
soz_loc_not = [squeeze(znew(1,1,:)),nanmean(squeeze([znew(2,2,:),znew(3,3,:)]),1)'];

figure
set(gcf,'position',[221 425 1220 452])
tiledlayout(1,2)
nexttile
stats = plot_paired_data(soz_lat_not',{'side of SOZ','side of non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);
title('Normalized intrinsic connectivity by SOZ laterality')

nexttile
stats = plot_paired_data(soz_loc_not',{'SOZ','non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);
title({'Normalized intrinsic connectivity by SOZ localization','(within same laterality)'})

print(gcf,[out_folder,'localize'],'-dpng');



end