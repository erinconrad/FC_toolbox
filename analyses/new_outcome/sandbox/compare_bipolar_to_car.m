function compare_bipolar_to_car

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;

%% Get variables of interest
pc_car = data.all_fc;
pc_bi = data.all_bipolar_fc;
coh_car = data.all_coh;
coh_bi = data.all_bipolar_coh;
bp_car = data.all_bp;
bp_bi = data.all_bipolar_bp;

%% plot a test couple
if 0
    figure
    nexttile
    turn_nans_gray(coh_car{100}(:,:,1))
    nexttile
    turn_nans_gray(coh_bi{100}(:,:,1))
end

%% Re-wrap
pc_car = cellfun(@wrap_or_unwrap_adjacency_fc_toolbox,...
    pc_car,'uniformoutput',false);
pc_bi = cellfun(@wrap_or_unwrap_adjacency_fc_toolbox,...
    pc_bi,'uniformoutput',false);
coh_car = cellfun(@wrap_or_unwrap_adjacency_fc_toolbox,...
    coh_car,'uniformoutput',false);
coh_bi = cellfun(@wrap_or_unwrap_adjacency_fc_toolbox,...
    coh_bi,'uniformoutput',false);

%% Get correlation between bipolar and car
r_pc = cellfun(@(x,y) ...
    corr(x,y,'rows','pairwise'),...
    pc_car,pc_bi);
r_coh = cellfun(@(x,y) ...
    (diag(corr(x,y,'rows','pairwise')))',...
    coh_car,coh_bi,'UniformOutput',false);
r_coh = cell2mat(r_coh);
r_bp = cellfun(@(x,y) ...
    (diag(corr(x,y,'rows','pairwise')))',...
    bp_car,bp_bi,'UniformOutput',false);
r_bp = cell2mat(r_bp);

plot(r_pc)

end