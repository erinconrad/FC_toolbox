function localize_epilepsy

%% Parameters
which_atlas = 'brainnetome';% %'aal';'aal_bernabei';
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


%% Get some basic measurements of connectivity
intrinsic = nan(nb,npts);
intrinsic_raw = nan(nb,npts);
extrinsic = nan(nb,npts);
for ip = 1:npts
    curr_pt = zbroad(:,:,ip);
    curr_pt_raw = broad_connectivity(:,:,ip);
    curr_intrinsic = curr_pt(logical(eye(size(curr_pt))));
    curr_intrinsic_raw = curr_pt_raw(logical(eye(size(curr_pt))));
    
    curr_pt_for_extrinsic = curr_pt;
    curr_pt_for_extrinsic(logical(eye(size(curr_pt)))) = nan;
    curr_extrinsic = nanmean(curr_pt_for_extrinsic,2);
    intrinsic(:,ip) = curr_intrinsic;
    intrinsic_raw(:,ip) = curr_intrinsic_raw;
    extrinsic(:,ip) = curr_extrinsic;
end

if 1
    figure
    set(gcf,'position',[100 100 1400 300])
    tiledlayout(1,5,'padding','tight','tilespacing','compact')
    
    nexttile
    turn_nans_gray(nanmean(broad_connectivity,3))
    %xlabel('Region')
    %ylabel('Region')
    title('Mean raw regional connectivity')
    yticks(1:nb)
    yticklabels(broad_regions)
    xticks(1:nb)
    xticklabels(broad_regions)
    set(gca,'fontsize',15)
    
    nexttile
    turn_nans_gray(nelecs_broad)
    xlabel('Patient')
    %ylabel('Region')
    title('# Elecs')
    yticklabels([])
    set(gca,'fontsize',15)
    
    nexttile
    turn_nans_gray(intrinsic)
    xlabel('Patient')
    %ylabel('Region')
    yticklabels([])
    title('Normalized intra-regional connectivity')
    set(gca,'fontsize',15)
    
    nexttile
    turn_nans_gray(extrinsic)
    xlabel('Patient')
    %ylabel('Region')
    yticklabels([])
    title('Normalized extra-regional connectivity')
    set(gca,'fontsize',15)
    
    nexttile
    turn_nans_gray(soz_broad)
    xlabel('Patient')
    %ylabel('Region')
    title('SOZ')
    yticklabels([])
    set(gca,'fontsize',15)
    
    
    
    print(gcf,[out_folder,'loc_methods_',which_atlas],'-dpng')
end

%% Get left/right connectivity
lr_conn = nan(npts,2);
left_idx = find(strcmp(broad_lat,'left'));
right_idx = find(strcmp(broad_lat,'right'));
soz_non_conn = nan(npts,2);
for ip = 1:npts
    % left conn
    left_conns = [];
    for ib = left_idx
        for jb = left_idx
            left_conns = [left_conns;(zbroad(ib,jb,ip))];
        end
    end
    lr_conn(ip,1) = nanmean(left_conns);
    
    % right conn
    right_conns = [];
    for ib = right_idx
        for jb = right_idx
            right_conns = [right_conns;(zbroad(ib,jb,ip))];
        end
    end
    lr_conn(ip,2) = nanmean(right_conns);
    
    curr_soz_broad = soz_broad(:,ip);
    if sum(curr_soz_broad) == 0 % skip if no SOZ (half the patients!)
        continue;
    end
    
    if ismember(find(curr_soz_broad),[1 3 5]) % left laterality
        soz_non_conn(ip,1) = lr_conn(ip,1);
        soz_non_conn(ip,2) = lr_conn(ip,2);
    elseif ismember(find(curr_soz_broad),[2 4 6]) % right laterality
        soz_non_conn(ip,2) = lr_conn(ip,1);
        soz_non_conn(ip,1) = lr_conn(ip,2);
    end
end


%% For each patient, compare the normalized regional connectivity of SOZ region to non-SOZ regions
soz_not = nan(npts,2,2); % first split for intrinsic vs extrinsic, second for soz vs not
%soz_not_non_normalized = nan(npts,2);
soz_lowest_n = nan(npts,2);
percentile = nan(npts,1);
for ip = 1:npts
    curr_soz_broad = soz_broad(:,ip); % get this patient's SOZ region
    
    if sum(curr_soz_broad) == 0 % skip if no SOZ (half the patients!)
        continue;
    end
    soz_not(ip,1,1) = (intrinsic(curr_soz_broad,ip)); % get normalized connectivity of SOZ
    soz_not(ip,1,2) = nanmean(intrinsic(~curr_soz_broad,ip)); % normalized connectivity of other regions, averaged across regions
    
    soz_not(ip,2,1) = (extrinsic(curr_soz_broad,ip)); % get normalized connectivity of SOZ
    soz_not(ip,2,2) = nanmean(extrinsic(~curr_soz_broad,ip)); % normalized connectivity of other regions, averaged across regions
    
    %{
    soz_not_non_normalized(ip,:) = [broad_connectivity(curr_soz_broad,ip),...
        nanmean(broad_connectivity(~curr_soz_broad,ip))]; %same but not normalized connectivity (not controlling for normal anatomical differences)
    %}
    
    %% get all zs for this patient
    %
    zcurr = intrinsic(:,ip);
    
    %% Get percentile (amongst non-nans) of SOZ
    zcurr_nans = isnan(zcurr);
    
    if sum(zcurr_nans) == length(zcurr), continue; end
    
    zcurr_non_nans = zcurr(~zcurr_nans);
    curr_soz_broad_non_nans = curr_soz_broad(~zcurr_nans);
    
    if ~any(curr_soz_broad_non_nans) == 1, continue; end
    
    [~,I] = sort(zcurr_non_nans);
    r = 1:length(zcurr_non_nans);
    r(I) = r;
    percentile(ip) = (r(curr_soz_broad_non_nans)+r(curr_soz_broad_non_nans)-1)/2/length(zcurr_non_nans);
    %}
    %% How often is SOZ lowest 
    [~,I] = min(zcurr_non_nans);
    if I == find(curr_soz_broad_non_nans)
        soz_lowest_n(ip,1) = 1;
    else
        soz_lowest_n(ip,1) = 0;
    end
    soz_lowest_n(ip,2) = length(zcurr_non_nans);
end

nan_rows = any(isnan(soz_lowest_n),2) | soz_lowest_n(:,2) == 0;
soz_lowest_n(nan_rows,:) = []; % note, I may be shooting myself in the foot by not removing those where soz is nan
percentile(nan_rows) = [];

%% Does soz have lowest connectivity more often than chance?
succ_and_p = [soz_lowest_n(:,1),1./soz_lowest_n(:,2)];
[pval,nchance] = general_binomial_test(succ_and_p);

%% For each localization (non-laterlized), compare connectivity between side with most electrodes and side with least
nbn = nb/2;
agree_conn_nelecs = nan(nbn,npts);
for ip = 1:npts
    for ib = 1:nbn
        curr_b = (ib-1)*2+1;
        
        if isnan(nelecs_broad(curr_b,ip)) || nelecs_broad(curr_b,ip) == nelecs_broad(curr_b+1,ip) || ...
                isnan(intrinsic_raw(curr_b,ip)) || isnan(intrinsic_raw(curr_b+1,ip)) || ...
                intrinsic_raw(curr_b,ip) == intrinsic_raw(curr_b+1,ip)
            continue;
        end
        more_elecs_left = nelecs_broad(curr_b,ip) > nelecs_broad(curr_b+1,ip);
        higher_conn_left = intrinsic_raw(curr_b,ip) > intrinsic_raw(curr_b+1,ip);
        
        if (more_elecs_left && higher_conn_left) || (~more_elecs_left && ~higher_conn_left)
            agree_conn_nelecs(ib,ip) = 1;
        else
            agree_conn_nelecs(ib,ip) = 0;
        end
    end
end

perc_agree = sum(agree_conn_nelecs == 1,2)./sum(~isnan(agree_conn_nelecs),2)*100;
new_ylabs = cell(3,1);
for i = 1:3
    new_ylabs{i} = sprintf('%s (%1.1f%%)',broad_non_lat{i},perc_agree(i));
end

figure
turn_nans_gray(agree_conn_nelecs)
yticks(1:3)
yticklabels(new_ylabs)
xticklabels([])
xlabel('Patient')
title('Does side with more electrodes have higher intrinsic connectivity?')

%% Localize epilepsy
figure
set(gcf,'position',[221 425 1220 372])
tiledlayout(1,2)

nexttile
stats = plot_paired_data(squeeze(soz_not(:,1,:))',{'SOZ','non-SOZ','non-SOZ'},'Normalized intrinsic connectivity','paired',plot_type);
title('Normalized intrinsic connectivity by SOZ status')

nexttile
stats = plot_paired_data(squeeze(soz_not(:,2,:))',{'SOZ','non-SOZ','non-SOZ'},'Normalized extrinsic connectivity','paired',plot_type);
title('Normalized extrinsic connectivity by SOZ status')

print(gcf,[out_folder,'loc_',which_atlas],'-dpng')

out.broad_connectivity = broad_connectivity;
out.soz_broad = soz_broad;
out.broad_regions = broad_regions;
save([out_folder,'broad_connectivity.mat'],'out');


%% Test if on an atlas-region level, the NS is different for SOZ than others
% Get node strengths (mean across edges for each node).
% I do mean instead of sum because nansum of all nans is 0, but nanmean of
% all nans is nan
if 0
ns = squeeze(nanmean(atlas,1));

% Convert ns to z-score, comparing each region to that of the same region in other patients
z = (ns - nanmean(ns,2))./nanstd(ns,[],2);

% For each patient, get the average z score of the SOZ electrodes and non-SOZ electrodes

z_soz_no = nan(npts,2);
for ip = 1:npts
    curr_z = z(:,ip);
    curr_sozs = sozs{ip};
    is_soz = ismember(atlas_nums,curr_sozs);
    z_soz = nanmean(curr_z(is_soz));
    z_not = nanmean(curr_z(~is_soz));
    z_soz_no(ip,:) = [z_soz z_not];
    
    
    if 0
        fprintf('\n%s SOZ:\n',pt_names{ip});
        names(is_soz)
        pause
        
    end
    
end

% Remove rows with any nans
any_nans = any(isnan(z_soz_no),2);
z_soz_no(any_nans,:) = [];

stats = plot_paired_data(z_soz_no',{'SOZ','Non-SOZ','non-SOZ'},'Node strength (z-score)','paired',plot_type);

fprintf(['\nThe mean edge weights of the SOZ regions (median %1.2f) was '...
    'lower than that of the non-SOZ regions (median %1.2f) (T+ = %1.2f, p = %1.3f).\n'],...
    stats.medians(1),stats.medians(2),stats.Tpos,stats.pval);

fprintf(['\nThe connectivity was higher in the non-SOZ for %d (%1.1f%%) of patients.\n'],...
    sum(diff(z_soz_no,1,2)>0),sum(diff(z_soz_no,1,2)>0)/size(z_soz_no,1)*100);
end


end