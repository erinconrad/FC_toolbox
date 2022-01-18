function lateralize_epilepsy

which_atlas = 'brainnetome';
plot_type = 'scatter';

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
nregions = length(names);
assert(nregions==size(atlas,1))

%% Get locs and lats for atlas names
[locs,lats,loc_nums] = lateralize_regions(names,which_atlas);
assert(isequal(lats,repmat({'L';'R'},123,1)))

%% Load soz lats
soz_out = load('out.mat');
soz_out = soz_out.out.circ_out;
soz_lats = soz_out.all_lats;
soz_pt_names = soz_out.names;
npts = length(pt_names);

% check names match
assert(isequal(pt_names,soz_pt_names))

%% Find regions in each patient that are bilaterally filled
all_bilateral = nan(nregions,npts);
for ip = 1:npts
    curr_atlas = atlas(:,:,ip);
    bilateral_region = zeros(nregions,1);
    for ir = 1:2:nregions-1
        
        % if both it and the next region (the contralateral one) have some
        % non-empty elements, keep it
        if sum(~isnan(curr_atlas(ir,:))) > 0 && ...
                sum(~isnan(curr_atlas(ir+1,:))) > 0
            bilateral_region(ir:ir+1) = [1;1];
        end
    end
    all_bilateral(:,ip) = bilateral_region;
end


%% Find patients with any bilateral
any_bilateral = (any(all_bilateral,1))';

%% Go through and check the regions make sense
% Look at a selection of patients and compare the regions that I think are
% bilateral to electrode labels for those patients
if 0
    %pt_names(any_bilateral)
    %ex = 'HUP100';
    %pt_idx = strcmp(pt_names,ex);
    bilat_names = cell(npts,1);
    
    for p = 1:npts
        bilat_names{p} = names(logical(all_bilateral(:,p)));
        %pt_names(p)
        %pause
    end
    
    all_lengths = cellfun(@length,bilat_names);
    max_length = max(all_lengths);
    for p = 1:npts
        bilat_names{p} = [bilat_names{p};repmat({''},max_length-all_lengths(p),1)];
    end
    
    C = cell(npts,length(bilat_names{1}));
    for i = 1:length(bilat_names)
        C(i,:) = bilat_names{i};
    end
    T = cell2table(C','VariableNames',pt_names);
    writetable(T,[out_folder,'bilateral_regions.csv'])
    
end

%% Find patients with unilateral epilepsy
unilateral_soz = strcmp(soz_lats,'left') | strcmp(soz_lats,'right');

%% Find patients for analysis: bilateral implants and unilateral epilepsy
include = any_bilateral & unilateral_soz;

%% Calculate intra-hemispheric FC for those to include
fc_lr = nan(npts,2);
% Loop over patients to include
for ip = 1:npts
    if include(ip) == 0, continue; end
    
    curr_atlas = atlas(:,:,ip);
    left = strcmp(lats,'L');
    right = strcmp(lats,'R');
    
    % Get the average intra-hemispheric edge weights
    left_fc = nanmean(curr_atlas(left,left),'all');
    right_fc = nanmean(curr_atlas(right,right),'all');
    
    fc_lr(ip,:) = [left_fc, right_fc];
end

%% Convert lr to soz-not soz
fc_soz_not = nan(npts,2);
for ip = 1:npts
    if include(ip) == 0, continue; end
    
    if strcmp(soz_lats{ip},'left')
        fc_soz_not(ip,:) = [fc_lr(ip,1) fc_lr(ip,2)];
    elseif strcmp(soz_lats{ip},'right')
        fc_soz_not(ip,:) = [fc_lr(ip,2) fc_lr(ip,1)];
    else
        error('what')
    end
end

%% remove nan rows
assert(isequal(any(isnan(fc_lr),2),any(isnan(fc_soz_not),2)))
any_nans = any(isnan(fc_lr),2);
fc_lr(any_nans,:) = [];
fc_soz_not(any_nans,:) = [];
nfinal = size(fc_lr,1);

%% Do some plots
p_soz = signrank(fc_soz_not(:,1),fc_soz_not(:,2));
p_lr = signrank(fc_lr(:,1),fc_lr(:,2));

figure
set(gcf,'position',[10 10 1000 400])
tiledlayout(1,2)

nexttile
plot_paired_data(fc_lr',{'left','right','right'},'Intra-hemispheric connectivity','paired',plot_type);
title('Left vs right connectivity')

nexttile
plot_paired_data(fc_soz_not',{'SOZ','non-SOZ','non-SOZ'},'Intra-hemispheric connectivity','paired',plot_type);
title('SOZ vs non-SOZ connectivity')

print(gcf,[out_folder,'lat'],'-dpng')

end