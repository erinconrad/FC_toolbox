function lateralize_epilepsy

which_atlas = 'aal_bernabei';%'brainnetome';%'aal_bernabei';% %
do_sw = 1;
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
nregions = length(names);
assert(nregions==size(atlas,1))

if contains(which_atlas,'aal_bernabei'), names = names'; end

%% Get locs and lats for atlas names
[locs,lats,loc_nums] = lateralize_regions(names,which_atlas);


%% Load soz lats
soz_out = load('out.mat');
soz_out = soz_out.out.circ_out;
soz_lats = soz_out.all_lats;
soz_pt_names = soz_out.names;
npts = length(pt_names);

% check names match
assert(isequal(pt_names,soz_pt_names))

%% Check some of these (looked good)
if 0
    table(soz_lats,soz_pt_names)
end

%% Find regions in each patient that are bilaterally filled
all_bilateral = nan(nregions,npts);

% get left regions
left_regions = strcmp(lats,'L');

for ip = 1:npts
    curr_atlas = atlas(:,:,ip);
    bilateral_region = zeros(nregions,1);
    
    % Loop over regions
    for ir = 1:nregions
        
        % skip if not a left region
        if left_regions(ir) == 0, continue; end
        
        % find the row of the corresponding right region
        curr_loc = locs{ir};
        right_region = strcmp(locs,curr_loc) & strcmp(lats,'R');
        
        % Skip if there is no corresponding right region
        if sum(right_region) == 0, continue; end
        
        % see if both the left and right region have non-empty elements in
        % the atlas
        if sum(~isnan(curr_atlas(ir,:))) > 0 && ...
                sum(~isnan(curr_atlas(right_region,:))) > 0
            bilateral_region(ir) = 1;
            bilateral_region(right_region) = 1;
        
        end
    end
    all_bilateral(:,ip) = bilateral_region;
end

%{
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
%}


%% Find patients with any bilateral
% Restrict to at least 4 bilateral so you can measure intrinsic
% connectivity on each side
%any_bilateral = (any(all_bilateral,1))';
any_bilateral = (sum(all_bilateral,1) >= 4)';

%% Go through and check the regions make sense
% Look at a selection of patients and compare the regions that I think are
% bilateral to electrode labels for those patients
if 0
    %pt_names(any_bilateral)
    %ex = 'HUP100';
    %pt_idx = strcmp(pt_names,ex);
    bilat_names = cell(npts,1);
    
    for p = 1:npts
        bilat_names{p} = (names(logical(all_bilateral(:,p))));
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
fprintf('\nThere are %d patients with both bilateral symmetric coverage and unilateral epilepsy\n',sum(include));

%% Calculate intra-hemispheric FC for those to include
% Could separately check intrinsic and extrinsic connectivity?
fc = cell(2,1);
fc{1} = nan(npts,2); % intrinsix
fc{2} = nan(npts,2); % extrinsic

fc_ws = cell(2,1);
fc_ws{1} = nan(npts,2,2);
fc_ws{2} = nan(npts,2,2);

% Loop over patients to include
for ip = 1:npts
    if include(ip) == 0, continue; end
    
    curr_atlas = atlas(:,:,ip);
    curr_atlas_ws = atlas_ws(:,:,:,ip);
    left = strcmp(lats,'L');
    right = strcmp(lats,'R');
    
    % Restrict to those regions with symmetric coverage
    curr_bilateral = logical(all_bilateral(:,ip));
    
    % Nice little code to visually check the restricted atlas for the
    % pateitn
    if 0
        
        % for visualization purposes, re-order by left and then right
        new_order = [find(left);find(right);find(~left&~right)];
        reordered_atlas = curr_atlas(new_order,new_order);
        reordered_names = names(new_order);
        reordered_bilateral = curr_bilateral(new_order);
        
        figure
        set(gcf,'position',[100 100 1100 350])
        tiledlayout(1,2)
        nexttile
        %imagesc(curr_atlas)
        turn_nans_gray(reordered_atlas)
        title(pt_names{ip})
        set(gca,'fontsize',15)
        xlabel('Parcels')
        ylabel('Parcels')
        %xticklabels([])
        %yticklabels([])
        
        nexttile
        turn_nans_gray(reordered_atlas(reordered_bilateral,reordered_bilateral))
        %imagesc(curr_atlas(curr_bilateral,curr_bilateral))
        xticks(1:sum(curr_bilateral))
        yticks(1:sum(curr_bilateral))
        xticklabels(reordered_names(reordered_bilateral))
        yticklabels(reordered_names(reordered_bilateral))
        set(gca,'fontsize',15)
        title('Symmetric parcellations')
        print(gcf,[out_folder,'lat_methods_',which_atlas],'-dpng')
        pause
        close(gcf)
    end
    
    % Get the average intrinsic edge weights of the symmetric regions
    %{
    This is to say, how connnected is each symmetric region to other
    symmetric regions on its same side? I think this is the most fair
    possible test because then I am only looking at the L sided regions
    that also have corresponding R sided regions.
    %}
    left_fc = nanmean(curr_atlas(left&curr_bilateral,left&curr_bilateral),'all');
    right_fc = nanmean(curr_atlas(right&curr_bilateral,right&curr_bilateral),'all');
    fc{1}(ip,:) = [left_fc, right_fc];
    
    % Get the average extrinsic edge weights of the symmetric regions
    %{
    This is to say, how connected is each symmetric region to all other
    regions?
    %}
    left_fc = nanmean(curr_atlas(left&curr_bilateral,~(left&curr_bilateral)),'all');
    right_fc = nanmean(curr_atlas(right&curr_bilateral,~(right&curr_bilateral)),'all');
    fc{2}(ip,:) = [left_fc, right_fc];
    
    % Intrinsic for ws
    left_fc_sw = nanmean(curr_atlas_ws(left&curr_bilateral,left&curr_bilateral,:),[1 2]);
    right_fc_sw = nanmean(curr_atlas_ws(right&curr_bilateral,right&curr_bilateral,:),[1 2]);
    fc_ws{1}(ip,:,1) = squeeze(left_fc_sw);
    fc_ws{1}(ip,:,2) = squeeze(right_fc_sw);
    
    % Extrinsic for ws
    left_fc_sw = nanmean(curr_atlas_ws(left&curr_bilateral,~(left&curr_bilateral),:),[1 2]);
    right_fc_sw = nanmean(curr_atlas_ws(right&curr_bilateral,~(right&curr_bilateral),:),[1 2]);
    fc_ws{2}(ip,:,1) = squeeze(left_fc_sw);
    fc_ws{2}(ip,:,2) = squeeze(right_fc_sw);
    
end

%% Convert lr to soz-not soz
fc_soz = cell(length(fc),1);
fc_soz{1} = nan(npts,2); % intrinsic
fc_soz{2} = nan(npts,2); % extrinsic

fc_soz_ws = cell(length(fc),1);
fc_soz_ws{1} = nan(npts,2,2);
fc_soz_ws{2} = nan(npts,2,2);
for i = 1:2
    fc_lr = fc{i};
    fc_ws_lr = fc_ws{i};
    for ip = 1:npts
        if include(ip) == 0, continue; end

        if strcmp(soz_lats{ip},'left')
            fc_soz{i}(ip,:) = [fc_lr(ip,1) fc_lr(ip,2)];
            fc_soz_ws{i}(ip,:,1) = squeeze(fc_ws_lr(ip,:,1));
            fc_soz_ws{i}(ip,:,2) = squeeze(fc_ws_lr(ip,:,2));
        elseif strcmp(soz_lats{ip},'right')
            fc_soz{i}(ip,:) = [fc_lr(ip,2) fc_lr(ip,1)];
            fc_soz_ws{i}(ip,:,1) = squeeze(fc_ws_lr(ip,:,2));
            fc_soz_ws{i}(ip,:,2) = squeeze(fc_ws_lr(ip,:,1));
        else
            error('what')
        end
    end
end

%% remove nan rows
for i = 1:2
    assert(isequal(any(isnan(fc{i}),2),any(isnan(fc_soz{i}),2)))
    any_nans = any(isnan(fc{i}),2);

    fc{i}(any_nans,:) = [];
    fc_soz{i}(any_nans,:) = [];

end

%% Also remove nans from ws
for i = 1:2
    % Think about this
end

if do_sw
    wake_fc_soz = squeeze(fc_soz_ws{1}(:,1,:));
    sleep_fc_soz = squeeze(fc_soz_ws{1}(:,2,:));
    
    figure
    set(gcf,'position',[10 10 1000 350])
    tiledlayout(1,2)
    nexttile
    stats = plot_paired_data(wake_fc_soz',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
    title('SOZ vs non-SOZ intrinsic connectivity: wake')
    
    nexttile
    stats = plot_paired_data(sleep_fc_soz',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
    title('SOZ vs non-SOZ intrinsic connectivity: sleep')
    
    
    
end

%% Do some plots
p_soz_in = signrank(fc_soz{1}(:,1),fc_soz{1}(:,2));
p_lr_in = signrank(fc{1}(:,1),fc{1}(:,2));
p_soz_ex = signrank(fc_soz{2}(:,1),fc_soz{2}(:,2));
p_lr_ex = signrank(fc{2}(:,1),fc{2}(:,2));

figure
set(gcf,'position',[10 10 1000 350])
tiledlayout(1,2)

nexttile
stats = plot_paired_data((fc{1})',{'left','right','right'},'Intrinsic connectivity','paired',plot_type);
title('Left vs right intrinsic connectivity')

nexttile
stats = plot_paired_data((fc_soz{1})',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
title('SOZ vs non-SOZ intrinsic connectivity')

%% text
fprintf(['\nThe intrinsic connectivity of the SOZ hemisphere (median %1.2f) was '...
    'lower than that of the non-SOZ hemisphere (median %1.2f) (T+ = %1.2f, p = %1.3f).\n'],...
    stats.medians(1),stats.medians(2),stats.Tpos,stats.pval);

fprintf(['\nThe connectivity was higher in the non-SOZ for %d (%1.1f%%) of patients.\n'],...
    sum(fc_soz{1}(:,1)-fc_soz{1}(:,2)<0),sum(fc_soz{1}(:,1)-fc_soz{1}(:,2)<0)/length(fc_soz{1})*100);


print(gcf,[out_folder,'lat_',which_atlas],'-dpng')


end