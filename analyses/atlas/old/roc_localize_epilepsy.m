function roc_localize_epilepsy

%% Parameters
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

%% Load intraregional connectivity matrix
out = load([out_folder,'broad_connectivity.mat']);
out = out.out;

bc = out.broad_connectivity;
soz = out.soz_broad;

% normalize across regions
z = (bc-nanmean(bc,2))./nanstd(bc,[],2);

% also normalize across patients
z = (z-nanmean(z,1))./nanstd(z,[],1);

% Vectorize
z = z(:);
soz = soz(:);

% remove nanrows
isnan_either = isnan(z) | isnan(soz);
z(isnan_either) = [];
soz(isnan_either) = [];

soz_label = cell(length(soz),1);
soz_label(soz) = {'soz'};
soz_label(~soz) = {'non-soz'};
pos_class = 'non-soz';
[X,Y,T,AUC,OPTROCPT] = perfcurve(soz_label,z,pos_class);
plot(X,Y)

end