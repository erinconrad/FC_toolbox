function localize_soz

%% Parameters
do_plots = 1;
which_atlas = 'aal_bernabei';%'aal_bernabei';%%'brainnetome';
plot_type = 'scatter';

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
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');

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


%% Does FC correlate with spike rate
z = (atlas-nanmean(atlas,3))./nanstd(atlas,[],3);
%ns = squeeze(nanmean(atlas,2));
ns = squeeze(nanmean(z,2));
if 0
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


% Make matrix identifying patient
which_pt = repmat(1:size(spikes,2),size(spikes,1),1);

% Vectorize spikes and node strength and which_pt
spike_vec = spikes(:);
ns_vec = ns(:);
pt_vec = which_pt(:);
soz_vec = bin_soz(:);

T = table(soz_vec,spike_vec,ns_vec,pt_vec);
T.pt_vec = nominal(T.pt_vec);
glme_raw = fitglme(T,'soz_vec ~ spike_vec + ns_vec + (1|pt_vec)',...
    'Distribution','Poisson','Link','log');

end