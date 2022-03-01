function fine_localization

%{
At individual region level, get node strength. Then normalize across
patients. Do SOZ regions have higher or lower connectivity than non SOZ
regions?

Edit this to:
- see if average intrinsic connectivity of Left different from right - no,
same
- see if average intrinsic connectivity of SOZ laterality different from
that of non-SOZ (don't normalize)
- then restrict to correct laterality.
- then exclude regions without many patients to normalize
- then normalize
- see if node strength different for SOZ vs non SOZ locations
%} 

%% Parameters
norm_pca = 0;
delete_nans = 1;
norm_with_normal =0;
do_plots = 0;
which_atlas = 'aal_bernabei';%%'brainnetome';
plot_type = 'scatter';
coverage_limit = 20;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];

bct_folder= locations.bct;
out_folder = [results_folder,'analysis/atlas/'];
plot_folder = [results_folder,'analysis/atlas/plots/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end
if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Get spikes
%{
sleep_data_folder = [scripts_folder,'analyses/sleep/data/'];
spikes_out = load([sleep_data_folder,'out.mat']);
spikes_out = spikes_out.out;
all_elecs_rates = spikes_out.bin_out.all_elecs_rates;
%}

%% Load atlas
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;

%% Load soz lats
%{
soz_out = load('out.mat');
soz_out = soz_out.out.circ_out;
soz_lats = soz_out.all_lats;
soz_locs = soz_out.all_locs;
%}
soz_locs = out.all_soz_locs;
soz_lats = out.all_soz_lats;

right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');

atlas = out.atlas;
atlas_norm = out.normal_atlas;
names = out.atlas_names;
pt_names = out.pt_names;
atlas_nums = out.atlas_nums;
nregions = length(names);
assert(nregions==size(atlas,1))
npts = size(atlas,3);
nelecs = out.n_elecs_all;
sozs = out.sozs;
bin_soz = (cell2mat(cellfun(@(x) ismember(atlas_nums',x),sozs,'uniformoutput',false)))';

%% Get spike rates in each parcel
spikes = out.spikes_atlas;
if 0
    figure
    turn_nans_gray(spikes)
    yticks(1:nregions)
    yticklabels(names)
end

%% Localize SOZ
mt = strcmp(soz_locs,'mesial temporal');
tn = strcmp(soz_locs,'temporal neocortical');
oc = strcmp(soz_locs,'other cortex');
temporal = contains(soz_locs,'temporal');
extra = strcmp(soz_locs,'other cortex') | strcmp(soz_locs,'diffuse') | strcmp(soz_locs,'multifocal');

%% Get locs and lats for atlas names
[locs,lats] = lateralize_regions(names,which_atlas);
left = strcmp(lats,'L');
right = strcmp(lats,'R');
lr_order = reorder_lr(locs,lats);

broad_locs = localize_regions(names,which_atlas);
broad_regions = {'left mesial temporal','right mesial temporal',...
    'left temporal neocortical','right temporal neocortical',...
    'left other cortex','right other cortex'};
nbroad = length(broad_regions);

%% Change all orders to left and right
left = left(lr_order);
right = right(lr_order);
atlas = atlas(lr_order,lr_order,:);
atlas_norm = atlas_norm(lr_order,lr_order,:);
names = names(lr_order);
atlas_nums = atlas_nums(lr_order);
bin_soz = bin_soz(lr_order,:);
spikes = spikes(lr_order,:);
broad_locs = broad_locs(lr_order);
locs = locs(lr_order);
lats = lats(lr_order);

%% Get left and right intrinsic connectivity
left_intrinsic = squeeze(nanmean(atlas(left,left,:),[1 2]));
right_intrinsic = squeeze(nanmean(atlas(right,right,:),[1 2]));
left_norm_intrinsic = squeeze(nanmean(atlas_norm(left,left,:),[1 2]));
right_norm_intrinsic = squeeze(nanmean(atlas_norm(right,right,:),[1 2]));

bl_intrinsic = [left_intrinsic right_intrinsic];

%% Sanity check (I expect no): Is left intrinsic connectivity diff from right?
% No
if do_plots
figure
stats = plot_paired_data((bl_intrinsic)',{'left','right','right'},'Intrinsic connectivity','paired',plot_type);
title('Left vs right intrinsic connectivity')
end

%% Get SOZ and non-SOZ laterality intrinsic connectivity
soz_non_intrinsic = nan(npts,2);
soz_non_intrinsic_norm = nan(npts,2);
for ip = 1:npts
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    
    if curr_soz_right == 1
        soz_non_intrinsic(ip,:) = [right_intrinsic(ip) left_intrinsic(ip)];
        soz_non_intrinsic_norm(ip,:) = [right_norm_intrinsic(ip) left_norm_intrinsic(ip)];
    elseif curr_soz_left == 1
        soz_non_intrinsic(ip,:) = [left_intrinsic(ip) right_intrinsic(ip)];
        soz_non_intrinsic_norm(ip,:) = [left_norm_intrinsic(ip) right_norm_intrinsic(ip)];
    end
end

%% I expect yes: Is SOZ intrinsic connectivity diff from non-SOZ?
% Yes
if do_plots
figure
tiledlayout(1,2)
nexttile
stats = plot_paired_data(soz_non_intrinsic',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
title('SOZ vs non-SOZ intrinsic connectivity')

nexttile
stats = plot_paired_data(soz_non_intrinsic_norm',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
title('SOZ vs non-SOZ intrinsic connectivity (using normal data)')
end

%% Confirmation of above test 1: for regions with symmetric coverage, is SOZ side intrinsic connectivity less?

%% First, build symmetric coverage atlas
symm_cov_atlas = nan(size(atlas));
all_bilateral = find_symmetric_coverage(atlas,lats,locs);
symm_soz_not = nan(npts,2);
for ip = 1:npts
    curr_bilateral = logical(all_bilateral(:,ip));
    if sum(curr_bilateral) == 0, continue; end
    curr_atlas = atlas(:,:,ip);
    curr_atlas(~curr_bilateral,:) = nan;
    curr_atlas(:,~curr_bilateral) = nan;
    symm_cov_atlas(:,:,ip) = curr_atlas;
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    if curr_soz_right == 1
        symm_soz_not(ip,:) = [nanmean(curr_atlas(right,right),'all') nanmean(curr_atlas(left,left),'all')];
    elseif curr_soz_left == 1
        symm_soz_not(ip,:) = [nanmean(curr_atlas(left,left),'all') nanmean(curr_atlas(right,right),'all')];
    end
end


%% Double check symm cov atlas is what I think


if strcmp(which_atlas,'brainnetome')
top_bilat = all_bilateral(1:123,:);
bot_bilat = all_bilateral(124:end,:);
assert(isequal(top_bilat,bot_bilat))
    
for ip = 1:npts
    curr = symm_cov_atlas(:,:,ip);
    for ir = 1:size(curr,1)/2
        if any(~isnan(curr(ir,:)))
            if ~any(~isnan(curr(ir+size(curr,1)/2,:)))
                error('what');
            end
        end
    end
end
end

if 0
    figure
    turn_nans_gray(all_bilateral)
end

% Is SOZ laterality different from non-SOZ laterality?
% Yes
if do_plots
figure
stats = plot_paired_data(symm_soz_not',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);
title({'SOZ vs non-SOZ intrinsic connectivity','(Symmetric coverage only)'})
end



%% For visualization, build a SOZ - non SOZ laterality ordered atlas
% Both for main and for symmetric coverage
soz_non_soz_ordered_atlas = nan(size(atlas));
soz_non_soz_ordered_atlas_symm_cov = nan(size(atlas));
for ip = 1:npts
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    curr_atlas = atlas(:,:,ip);
    curr_symm = symm_cov_atlas(:,:,ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    
    if curr_soz_right == 1
        soz_order = [find(right);find(left)];
    elseif curr_soz_left == 1
        soz_order = [find(left);find(right)];
    end
    
    soz_non_soz_ordered_atlas(1:length(soz_order),1:length(soz_order),ip) = curr_atlas(soz_order,soz_order);
    soz_non_soz_ordered_atlas_symm_cov(1:length(soz_order),1:length(soz_order),ip) = curr_symm(soz_order,soz_order);
    
end

names_minus_lat = cellfun(@(x) x(1:end-9),names,'uniformoutput',false);

%% Show my SOZ-non SOZ ordered matrices
if do_plots
    figure
    set(gcf,'position',[10 10 1000 500])
    tiledlayout(1,2,'tilespacing','compact','padding','tight')
    
    nexttile
    turn_nans_gray(nanmean(soz_non_soz_ordered_atlas,3))
    yticks(1:3:length(names))
    yticklabels(names_minus_lat(1:3:length(names)))
    set(gca,'fontsize',15)
    xticklabels([])
    title('All regions')
    
    nexttile
    turn_nans_gray(nanmean(soz_non_soz_ordered_atlas_symm_cov,3))
    xticklabels([])
    yticklabels([])
    title('Symmetric coverage only')
    set(gca,'fontsize',15)
    
end

%% Show number of patients contributing to each area
cov_map = nan(size(atlas,1),size(atlas,2));
symm_cov_map = nan(size(atlas,1),size(atlas,2));

for ir = 1:size(atlas,1)
    for jr = 1:size(atlas,2)
        cov_map(ir,jr) = sum(~isnan(soz_non_soz_ordered_atlas(ir,jr,:)));
        symm_cov_map(ir,jr) = sum(~isnan(soz_non_soz_ordered_atlas_symm_cov(ir,jr,:)));
    end
end

if do_plots
    figure
    set(gcf,'position',[10 10 1000 500])
    tiledlayout(1,2,'tilespacing','compact','padding','tight')
    
    nexttile
    turn_nans_gray((cov_map))
    yticks(1:5:length(names))
    yticklabels(names_minus_lat(1:5:length(names)))
    set(gca,'fontsize',15)
    xticklabels([])
    title('All regions - coverage')
    colorbar
    
    nexttile
    turn_nans_gray(symm_cov_map)
    xticklabels([])
    yticklabels([])
    title('Symmetric coverage only')
    set(gca,'fontsize',15)
    colorbar
end

%{
%% Alt symm coverage test
% I confirmed this gives the same result as the original  symm coverage
% test above
midpoint = size(atlas,1)/2;
lu_square = soz_non_soz_ordered_atlas_symm_cov(1:midpoint,1:midpoint,:);
rl_square = soz_non_soz_ordered_atlas_symm_cov(midpoint+1:end,midpoint+1:end,:);
lu_mean = squeeze(nanmean(lu_square,[1 2]));
rl_mean = squeeze(nanmean(rl_square,[1 2]));

if 0
    figure
    stats = plot_paired_data(([lu_mean rl_mean])',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity','paired',plot_type);

end
%}

%% Accuracy at lateralizing epilepsy
% Get average spike rate by laterality
spikes_lr = nan(npts,2);
spikes_lr(:,1) = nanmean(spikes(left,:));
spikes_lr(:,2) = nanmean(spikes(right,:));
spikes_higher_left = diff(spikes_lr,[],2) < 0;

% Get connectivity by left and right
fc_lr = nan(npts,2);
fc_lr(:,1) = nanmean(atlas(left,left,:),[1 2]);
fc_lr(:,2) = nanmean(atlas(right,right,:),[1 2]);
fc_lower_left = diff(fc_lr,[],2) > 0;

% Remove patients without epilepsy lateralization or any nans
no_lat = ~left_lat & ~right_lat;
any_nans = any(isnan(spikes_lr),2) | any(isnan(fc_lr),2);
to_remove = no_lat | any_nans;

left_rm = left_lat(~to_remove);
spikes_higher_left = spikes_higher_left(~to_remove);
fc_lower_left = fc_lower_left(~to_remove);

% Confusion matrix for each
spikes_conf = confusion_matrix(spikes_higher_left,left_rm,0);
fc_conf = confusion_matrix(fc_lower_left,left_rm,0);



%% Can I distinguish between bilateral and unilateral epilepsy?
bilat = strcmp(soz_lats,'bilateral');
unilat = strcmp(soz_lats,'left') | strcmp(soz_lats,'right');
temporal = contains(soz_locs,'temporal');

% Get transitivity of all patients
all_T = nan(npts,1);
for ip = 1:npts
    curr = symm_cov_atlas(:,:,ip);
    curr(isnan(curr)) = nanmean(curr,'all');
    all_T(ip)=transitivity_wu(curr);
end

% avg connectivity
all_avg = nan(npts,1);
all_intrinsic_avg = nan(npts,1);
for ip = 1:npts
    curr = symm_cov_atlas(:,:,ip);
    all_avg(ip) = nanmean(curr,'all');
    all_intrinsic_avg(ip) = nanmean([nanmean(curr(left,left),'all'),...
        nanmean(curr(right,right),'all')]);
    
end

% modularity
all_Q = nan(npts,1);
for ip = 1:npts
    curr = symm_cov_atlas(:,:,ip);
    curr(isnan(curr)) = nanmean(curr,'all');
    if any(isnan(curr),'all'), continue; end
    [~,all_Q(ip)]=modularity_und(curr,1);
end

% clust coeff
all_C = nan(npts,1);
for ip = 1:npts
    curr = symm_cov_atlas(:,:,ip);
    curr(isnan(curr)) = nanmean(curr,'all');
    if any(isnan(curr),'all'), continue; end
    C=clustering_coef_wu(curr);
    all_C(ip) = nanmean(C);
end


if 0
    figure
    tiledlayout(1,2)
    
    nexttile
    turn_nans_gray(nanmean(atlas(:,:,bilat),3))
    
    nexttile
    turn_nans_gray(nanmean(soz_non_soz_ordered_atlas(:,:,~bilat),3))
end

% average INTRINSIC connectivity (L-L, R-R)
% my rationale for doing this as opposed to overall avg connectivity is
% that bilateral patients should have broader coverage and potentially
% lower connectivity)

% show transitivity of bilat vs unilat
if 0
    figure
    plot(1+randn(sum(bilat),1)*0.05,all_T(bilat),'o','linewidth',2)
    hold on
    plot(2+randn(sum(unilat),1)*0.05,all_T(unilat),'o','linewidth',2)
end

%% prep stuff for lateralityfig
lat_fig.atlas = atlas;
lat_fig.soz_atlas = soz_non_soz_ordered_atlas;
lat_fig.symm_soz_atlas = soz_non_soz_ordered_atlas_symm_cov;
lat_fig.cov_map = cov_map;
lat_fig.symm_cov = symm_cov_map;
lat_fig.symm_soz_not = symm_soz_not;
lat_fig.bl_intrinsic = bl_intrinsic;
lat_fig.soz_non_intrinsic = soz_non_intrinsic;
lat_fig.spikes_conf = spikes_conf;
lat_fig.fc_conf = fc_conf;
lat_fig.names_minus_lat = names_minus_lat;
lat_fig.names = names;
lat_fig.Q = all_Q;
lat_fig.bilat = bilat;
lat_fig.unilat = unilat;
lat_fig.plot_folder = plot_folder;

laterality_fig(lat_fig)

%% Get average connectivity of broad regions
broad_conn = nan(nbroad,npts);
for ip = 1:npts
    curr_atlas = atlas(:,:,ip);
    for ib = 1:nbroad
        broad_conn(ib,ip) = squeeze(nanmean(curr_atlas(strcmp(broad_locs,broad_regions{ib}),:),'all'));
    end
end

% How does raw connectivity differ regionally?
if 0
    turn_nans_gray(broad_conn)
    yticks(1:nbroad)
    yticklabels(broad_regions)
end

%% For future analyses, remove regions with low coverage
% find regions with low coverage
if norm_with_normal
    any_coverage = squeeze(any(~isnan(atlas_norm),2));
else
    any_coverage = squeeze(any(~isnan(atlas),2));  
end

low_coverage = sum(any_coverage,2) < coverage_limit;
atlas(low_coverage,:,:) = nan;
atlas_norm(low_coverage,:,:) = nan;
atlas_norm(:,low_coverage,:) = nan;
atlas(:,low_coverage,:) = nan;
spikes(low_coverage) = nan;

%% Normalize edges
if norm_with_normal
    z = (atlas - nanmean(atlas_norm,3))./nanstd(atlas_norm,[],3);
else
    z = (atlas-nanmean(atlas,3))./nanstd(atlas,[],3);
end

%z = abs(z);

%% Confirmation of above test 2: is normalized SOZ intrinsic connectivity diff from non-SOZ?
% Get left and right intrinsic connectivity
left_intrinsic_z = squeeze(nanmean(z(left,left,:),[1 2]));
right_intrinsic_z = squeeze(nanmean(z(right,right,:),[1 2]));

% Sanity check (I expect no): Is left intrinsic connectivity diff from right?
% No
if do_plots
figure
stats = plot_paired_data(([left_intrinsic_z right_intrinsic_z])',{'left','right','right'},'Intrinsic connectivity','paired',plot_type);
title('Left vs right intrinsic connectivity (normalized)')
end

% Get SOZ and non-SOZ laterality intrinsic connectivity
soz_non_intrinsic_norm = nan(npts,2);
for ip = 1:npts
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    
    if curr_soz_right == 1
        soz_non_intrinsic_norm(ip,:) = [right_intrinsic_z(ip) left_intrinsic_z(ip)];
    elseif curr_soz_left == 1
        soz_non_intrinsic_norm(ip,:) = [left_intrinsic_z(ip) right_intrinsic_z(ip)];
    end
end

%  is normalized SOZ intrinsic connectivity diff from non-SOZ?
% yes
if do_plots
figure
stats = plot_paired_data(soz_non_intrinsic_norm',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity (norm)','paired',plot_type);
title('SOZ vs non-SOZ side intrinsic connectivity (normalized)')
end

%% Restrict to SOZ laterality, get SOZ vs non SOZ normalized connectivity
soz_non_total = nan(npts,2);
soz_non_in = nan(npts,2);
for ip = 1:npts
    
    % Get SOZ laterality
    curr_soz_right = right_lat(ip);
    curr_soz_left = left_lat(ip);
    
    if ~curr_soz_right && ~curr_soz_left, continue; end
    
    % Get current normalized atlas
    curr_z = z(:,:,ip);
    
    % Make non-soz side all nans
    if curr_soz_right
        %curr_z(left,:) = nan; curr_z(:,left) = nan;
    elseif curr_soz_left
        %curr_z(right,:) = nan; curr_z(:,right) = nan;
    end
    
    % get SOZ regions
    curr_soz_region = bin_soz(:,ip);
    
    % compare mean SOZ total connectivity to non
    soz_non_total(ip,1) = nanmean(curr_z(curr_soz_region,:),'all');
    soz_non_total(ip,2) = nanmean(curr_z(~curr_soz_region,:),'all');
    
    % compare mean SOZ intrinsic connectivity to non
    soz_non_in(ip,1) = nanmean(curr_z(curr_soz_region,curr_soz_region),'all');
    soz_non_in(ip,2) = nanmean(curr_z(~curr_soz_region,~curr_soz_region),'all');
    
end

%% Is SOZ different from same laterality non-SOZ?
% NO
if do_plots
figure
stats = plot_paired_data(soz_non_total',{'SOZ','non-SOZ','non-SOZ'},'Intrinsic connectivity (norm)','paired',plot_type);
title('SOZ vs non-SOZ intrinsic connectivity (normalized), same side')
end


%% Does NS correlate with spike rate?
% Yes! Higher average connectivity (raw or normalized) -> higher spike rate
ns_norm = squeeze(nanmean(z,2));
ns = squeeze(nanmean(atlas,2));

if 0
    figure
    for ip = 1:npts
        plot(ns(:,ip),spikes(:,ip),'wo')
        hold on
        text(ns(:,ip),spikes(:,ip),names,'horizontalalignment','center')
        xlabel('Mean connectivity')
        ylabel('Spike rate')
        pause
        hold off
    end
end

if do_plots
all_corrs = nan(npts,1);
figure
for ip = 1:npts
    all_corrs(ip) = corr(ns(:,ip),spikes(:,ip),'type','spearman','rows','pairwise');
    plot(ip,all_corrs(ip),'ko','linewidth',2)
    hold on
    
end
plot(xlim,[0 0],'k--','linewidth',2)
xlabel('Patient')
ylabel('Correlation between spike rate and connectivity')
end

%% Logistic regression mixed effects model: controlling for spike rate, is connectivity lower in SOZ?
% Yes, but only when using brainnetome atlas

% Make matrix identifying patient
which_pt = repmat(1:size(spikes,2),size(spikes,1),1);

% Vectorize spikes and node strength and which_pt
spike_vec = spikes(:);
ns_vec = ns(:);
ns_norm_vec = ns_norm(:);
pt_vec = which_pt(:);
soz_vec = bin_soz(:);

T = table(soz_vec,spike_vec,ns_vec,ns_norm_vec,pt_vec);

%% machine learning - can I improve localization of SOZ if I use FC rather than spikes alone? 
% NO
if 0
all_auc = nan(1e3,1);
for ib = 1:1e3
    
    all_auc(ib) = predict_soz_fc(T);
end
end

%% Logistic regression - If I take intrinsic connectivity and spikes, can I lateralize epilepsy


%% LR model

%{
T.pt_vec = nominal(T.pt_vec);
glme_raw = fitglme(T,'soz_vec ~ spike_vec + ns_vec + (1|pt_vec)',...
    'Distribution','Poisson','Link','log');

glme_norm = fitglme(T,'soz_vec ~ spike_vec + ns_norm_vec + (1|pt_vec)',...
    'Distribution','Poisson','Link','log');

glme_spikes = fitglme(T,'soz_vec ~ spike_vec + (1|pt_vec)',...
    'Distribution','Poisson','Link','log');

glme_ns_raw = fitglme(T,'soz_vec ~ ns_vec + (1|pt_vec)',...
    'Distribution','Poisson','Link','log');

glme_ns_norm = fitglme(T,'soz_vec ~ ns_norm_vec + (1|pt_vec)',...
    'Distribution','Poisson','Link','log');

%}
%{

%% Use normalized connectivity! (Because the biggest signal appears to be this non-normalized anatomical effect)
if norm_pca == 1
    atlas = z;
elseif norm_pca == 2
    atlas = (atlas - nanmean(atlas,[1 2]))./nanstd(atlas,[],[1 2]);
end



%% Remove nans
% Remove all nan observations

all_nan_rows = sum(isnan(nanmean(atlas,3)),1) == size(atlas,2);
atlas_rm = atlas;
names_rm = names;
left_rm = left;
right_rm = right;

if delete_nans

    atlas_rm(all_nan_rows,:,:) = [];
    atlas_rm(:,all_nan_rows,:) = [];
    names_rm(all_nan_rows) = [];
    left_rm(all_nan_rows) = [];
    right_rm(all_nan_rows) = [];

end

% Remove all nan patients
all_nan_pts = sum(isnan(atlas),[1 2]) == size(atlas,1)*size(atlas,2);
atlas_rm(:,:,all_nan_pts) = [];
right_non_nan = right_lat(~all_nan_pts);
left_non_nan = left_lat(~all_nan_pts);
tn = tn(~all_nan_pts);
mt = mt(~all_nan_pts);
oc = oc(~all_nan_pts);
temporal = temporal(~all_nan_pts);
extra = extra(~all_nan_pts);


% wrap
wrapped_atlas = (wrap_or_unwrap_adjacency_fc_toolbox(atlas_rm))';
nrows = size(wrapped_atlas,1);

%% Normalize
% If I do PCA without prenormalizing the data, the first component seems to
% be highly variable ACROSS patients but not too variable within patients.
% I think I should normalize within a patient.
%wrapped_atlas = (wrapped_atlas-nanmean(wrapped_atlas,[1 2]))./nanstd(wrapped_atlas,[],[1 2]);


%% Do NNMF
% Fill in nans with mean for that location
nnmf_atlas = atlas_rm;

for ip = 1:size(atlas_rm,3)
    %nnmf_atlas(logical(eye(size(nnmf_atlas(:,:,ip)))),ip) = 0;
    for i = 1:size(atlas_rm,1)
        for j = 1:size(atlas_rm,2)
            if i == j, nnmf_atlas(i,j,ip) = 0; continue; end
            if isnan(nnmf_atlas(i,j,ip))
                nnmf_atlas(i,j,ip) = nanmean(atlas_rm(i,j,:));
                
                % if no examples, make it the mean of that row for that
                % patient
                if isnan(nnmf_atlas(i,j,ip))
                    nnmf_atlas(i,j,ip) = nanmean(squeeze(nnmf_atlas(:,j,:)),'all');
                end
                
                if isnan(nnmf_atlas(i,j,ip)), error('why'); end
                    
            end
        end
    end
end

if 0
    figure
    turn_nans_gray(nanmean(nnmf_atlas,3))
end

% wrap it
nnmf_wrap = (wrap_or_unwrap_adjacency_fc_toolbox(nnmf_atlas))';

% NNMF
[W,H] = nnmf(nnmf_wrap,2);

if 1
    figure
    set(gcf,'position',[10 10 1400 500])
    tiledlayout(1,length(things))
    for i = 1:2
        nexttile
        uw = wrap_or_unwrap_adjacency_fc_toolbox(H(i,:)');
        
        uw_new = uw;
        names_new = names_rm;
        
        if ~delete_nans
            names_new = names_new(~all_nan_rows);
            uw_new = uw_new(~all_nan_rows,~all_nan_rows);
        end
        
        turn_nans_gray(uw_new)
        
        yticks(1:length(names_new))
        yticklabels(names_new)
        
        xticks(1:length(names_new))
        xticklabels(names_new)
        
    end
end

%% PCA
% PCA using ALS algorithm for missing data
[coeff,score,latent,~,explained,mu1] = pca(wrapped_atlas,'algorithm','als','Centered',false);


%% Put output into table for left right
non_lateralized = ~left_non_nan & ~right_non_nan;
tleft = left_non_nan;
tright = right_non_nan;
tscore = score;
tleft(non_lateralized) = [];
tright(non_lateralized) = [];
tscore(non_lateralized,:) = [];
T = array2table(tscore);
T = addvars(T,tleft,'Before','tscore1');

%% Show some diagnostics
if 0
    
    figure
    tiledlayout(1,2)
    nexttile
    turn_nans_gray(nanmean(atlas_rm,3))
    
    % rebuild data
    rebuild = score*coeff';
    rebuild_w = wrap_or_unwrap_adjacency_fc_toolbox(rebuild');
    nexttile
    turn_nans_gray(nanmean(rebuild_w,3))
end

if 0
    figure
    tiledlayout(1,2)
    nexttile
    turn_nans_gray(wrapped_atlas)
    
    
    nexttile
    rebuild = score*coeff';
    new_rebuild = rebuild;
    new_rebuild(isnan(wrapped_atlas)) = nan;
    turn_nans_gray(new_rebuild)
    caxis([0 1])
end


if 1
    figure
    set(gcf,'position',[10 10 1400 500])
    things = [1:4];
    tiledlayout(1,length(things))
    for i = things
        nexttile
        uw = wrap_or_unwrap_adjacency_fc_toolbox(coeff(:,i));
        
        uw_new = uw;
        names_new = names_rm;
        
        if ~delete_nans
            names_new = names_new(~all_nan_rows);
            uw_new = uw_new(~all_nan_rows,~all_nan_rows);
        end
        
        turn_nans_gray(uw_new)
        
        yticks(1:length(names_new))
        yticklabels(names_new)
        
        xticks(1:length(names_new))
        xticklabels(names_new)
        
    end
end

idx = kmeans(score(:,1:4),2);
if 0
    figure
    plot(score(idx==1,1),score(idx==1,2),'o')
    hold on
    plot(score(idx==2,1),score(idx==2,2),'o')
end

if 0
    figure
    plot(randn(sum(left_non_nan),1)*0.05+1,idx(left_non_nan),'o')
    hold on
    plot(randn(sum(right_non_nan),1)*0.05+2,idx(right_non_nan),'o')
end

if 0
    figure
    plot(randn(sum(mt),1)*0.05+1,idx(mt),'o')
    hold on
    plot(randn(sum(tn),1)*0.05+2,idx(tn),'o')
    plot(randn(sum(oc),1)*0.05+3,idx(oc),'o')
end

if 0
    which_scores = [1 2];
    figure
    plot(score(left_non_nan,which_scores(1)),score(left_non_nan,which_scores(2)),'o','linewidth',2)
    hold on
    plot(score(right_non_nan,which_scores(1)),score(right_non_nan,which_scores(2)),'o','linewidth',2)
    xlabel(sprintf('Score %d',which_scores(1)))
    ylabel(sprintf('Score %d',which_scores(2)))
    legend({'Left epilepsy','Right epilepsy'})
end

if 0
    which_scores = [1 4];
    figure
    plot(score(temporal,which_scores(1)),score(temporal,which_scores(2)),'o','linewidth',2)
    hold on
    plot(score(extra,which_scores(1)),score(extra,which_scores(2)),'o','linewidth',2)
    xlabel(sprintf('Score %d',which_scores(1)))
    ylabel(sprintf('Score %d',which_scores(2)))
    legend({'Temporal','Extra-temporal'})
end

if 0
    figure
    which_scores = [1 2];
    plot(score(mt,which_scores(1)),score(mt,which_scores(2)),'o','linewidth',2)
    hold on
    plot(score(tn,which_scores(1)),score(tn,which_scores(2)),'o','linewidth',2)
    plot(score(oc,which_scores(1)),score(oc,which_scores(2)),'o','linewidth',2)
    xlabel(sprintf('Score %d',which_scores(1)))
    ylabel(sprintf('Score %d',which_scores(2)))
    legend({'Mesial temporal','temporal neocortical','other cortex'})
end

%}

end