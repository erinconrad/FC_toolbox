function localize_epilepsy

%% Parameters
do_sw = 1;
which_atlas = 'brainnetome';%'aal_bernabei';%% %'aal';'aal_bernabei';
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
if do_sw
    out = load([atlas_folder,which_atlas,'_ws.mat']);
else
    out = load([atlas_folder,which_atlas,'.mat']);
end

out = out.out;

atlas = out.atlas;
atlas_ws = out.atlas_ws;
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
broad_connectivity_ws = nan(nb,nb,2,npts); % wake, sleep
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
            
            conn_ws = squeeze(nanmean(atlas_ws(i_regions,j_regions,:,ip),[1 2]));
            broad_connectivity_ws(ib,jb,:,ip) = conn_ws;

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
zbroad_ws = (broad_connectivity_ws - nanmean(broad_connectivity_ws,4))./nanstd(broad_connectivity_ws,[],4);


%% Reorder connectivity matrix by SOZ - non SOZ ipsi 1 - non SOZ ipsi 2 - contralateral for each
znew = reorder_locs_soz(zbroad,soz_broad,0);
if do_sw
    znew_ws = reorder_locs_soz(zbroad_ws,soz_broad,1);
end
if 0
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

% WRITE THIS CODE
soz_lat_not_wake = [nanmean(squeeze([znew_ws(1,1,1,:),znew_ws(2,2,1,:),znew_ws(3,3,1,:)]),1)',...
    nanmean(squeeze([znew_ws(4,4,1,:),znew_ws(5,5,1,:),znew_ws(6,6,1,:)]),1)'];
soz_loc_not_wake = [squeeze(znew_ws(1,1,1,:)),nanmean(squeeze([znew_ws(2,2,1,:),znew_ws(3,3,1,:)]),1)'];

soz_lat_not_sleep = [nanmean(squeeze([znew_ws(1,1,2,:),znew_ws(2,2,2,:),znew_ws(3,3,2,:)]),1)',...
    nanmean(squeeze([znew_ws(4,4,2,:),znew_ws(5,5,2,:),znew_ws(6,6,2,:)]),1)'];
soz_loc_not_sleep = [squeeze(znew_ws(1,1,2,:)),nanmean(squeeze([znew_ws(2,2,2,:),znew_ws(3,3,2,:)]),1)'];


%% Show sleep and wake
if do_sw
    figure
    set(gcf,'position',[221 425 1220 800])
    tiledlayout(2,2)
    nexttile
    stats = plot_paired_data(soz_lat_not_wake',{'side of SOZ','side of non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);
    title('Normalized intrinsic connectivity by SOZ laterality - wake')

    nexttile
    stats = plot_paired_data(soz_loc_not_wake',{'SOZ','non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);
    title({'Normalized intrinsic connectivity by SOZ localization - wake','(within same laterality)'})
    
    nexttile
    stats = plot_paired_data(soz_lat_not_sleep',{'side of SOZ','side of non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);
    title('Normalized intrinsic connectivity by SOZ laterality - sleep')

    nexttile
    stats = plot_paired_data(soz_loc_not_sleep',{'SOZ','non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);
    title({'Normalized intrinsic connectivity by SOZ localization - sleep','(within same laterality)'})
    
end

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