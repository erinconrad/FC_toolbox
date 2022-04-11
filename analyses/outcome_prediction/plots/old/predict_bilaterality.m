function predict_bilaterality

%% Parameters
do_plots = 1;
which_atlas = 'aal_bernabei';%'brainnetome';%%%
plot_type = 'scatter';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];

bct_folder= locations.bct;
out_folder = [results_folder,'analysis/atlas/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end
if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));


%% Load atlas and get region names and spikes
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;
atlas = out.atlas;
names = out.atlas_names;
spikes = out.spikes_atlas;
nregions = length(names);
assert(nregions==size(atlas,1))
npts = size(atlas,3);
sozs = out.sozs;
atlas_nums = out.atlas_nums;
bin_soz = (cell2mat(cellfun(@(x) ismember(atlas_nums',x),sozs,'uniformoutput',false)))';

%% Load soz lats
soz_lats = out.all_soz_lats;
bilat = cellfun(@(x) contains(x,'bilateral') || contains(x,'diffuse'),soz_lats);
unilat = cellfun(@(x) contains(x,'left') || contains(x,'right'),soz_lats);

%% Define names corresponding to mesial temporal
switch which_atlas
    case 'aal_bernabei'
        mt_names = {'Hippocampus','Amygdala'};
    case 'brainnetome'
        mt_names = {'Amyg','Hipp'};
end

mt = contains(names,mt_names);

%% Get locs and lats for atlas names
[locs,lats] = lateralize_regions(names,which_atlas);
left = strcmp(lats,'L');
right = strcmp(lats,'R');
neither_lat = ~left & ~right;

% confirm atlas has as many on right as left
assert(sum(left)==sum(right));

%% Re-order atlas to be left then right then neither
lr_order = reorder_lr(locs,lats);

left = left(lr_order);
right = right(lr_order);
neither_lat = neither_lat(lr_order);
atlas = atlas(lr_order,lr_order,:);
names = names(lr_order);
locs = locs(lr_order);
lats = lats(lr_order);
spikes = spikes(lr_order,:);
bin_soz = bin_soz(lr_order,:);
atlas_nums = atlas_nums(lr_order);
mt = mt(lr_order);

%% index of contralateral
contra_index = nan(size(left));
contra_index(1:sum(left)) = ([1:sum(left)])'+sum(left);
contra_index(sum(left)+1:sum(left)*2) = ([1:sum(left)])';

% make sure first half locs are same as last half locs
assert(isequal(locs(1:sum(left)),locs(sum(left)+1:sum(left)*2)))
assert(isequal(locs(contra_index(1:sum(left))),locs(contra_index(sum(left)+1:sum(left)*2))))

%% Build atlas ONLY containing mesial temporal regions
mt_atlas = atlas;
mt_atlas(~mt,:,:) = nan; mt_atlas(:,~mt,:) = nan;
mt_spikes = spikes;
mt_spikes(~mt,:) = nan;
mt_names = names;
mt_names(~mt) = {''};

%% Show stuff
if 0
    figure; tiledlayout(1,2);
    nexttile
    turn_nans_gray(nanmean(mt_atlas,3))
    xticks(1:length(mt_names))
    yticks(1:length(mt_names))
    xticklabels(mt_names)
    yticklabels(mt_names)
    
    nexttile
    turn_nans_gray(mt_spikes)
    yticks(1:length(mt_names))
    yticklabels(mt_names)
    xlabel('Patient')
end

%% Get average L-R connectivity for each patient
lr_conn = nan(npts,1);
lr_spikes = nan(npts,2);
for ip = 1:npts
    assert(nanmean(mt_atlas(left,right,ip),'all') == nanmean(mt_atlas(right,left,ip),'all') ...
        || isnan(nanmean(mt_atlas(left,right,ip),'all')))
    lr_conn(ip) = nanmean(mt_atlas(left,right,ip),'all');
    lr_spikes(ip,:) = [nanmean(mt_spikes(left,ip)) nanmean(mt_spikes(right,ip))];
end

%% Compare lr conn between unilat and bilateral
if 0
    figure
    plot(1+randn(sum(unilat),1)*0.05,lr_conn(unilat),'o','linewidth',2)
    hold on
    plot(2+randn(sum(bilat),1)*0.05,lr_conn(bilat),'o','linewidth',2)
    xticks([1 2])
    xticklabels({'Unilateral','bilateral'})
    ylabel('Left mesial temporal - right mesial temporal connectivity')
    p = ranksum(lr_conn(bilat),lr_conn(unilat));
    title(sprintf('%s',get_p_text(p)))
    set(gca,'fontsize',15)
end

%% Compare spikes left and right unilat and bilat
abs_diff = abs((lr_spikes(:,1)-lr_spikes(:,2))./nanmean(lr_spikes,2));
if 0
    figure
    plot(1+randn(sum(unilat),1)*0.05,abs_diff(unilat),'o','linewidth',2)
    hold on
    plot(2+randn(sum(bilat),1)*0.05,abs_diff(bilat),'o','linewidth',2)
    xticks([1 2])
    xticklabels({'Unilateral','bilateral'})
    ylabel('Unsigned relative difference in spikes unilateral vs bilateral')
    p = ranksum(abs_diff(bilat),abs_diff(unilat));
    title(sprintf('%s',get_p_text(p)))
    set(gca,'fontsize',15)
end


end