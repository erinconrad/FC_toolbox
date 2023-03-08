function outcome_plots(T,features)

%% Parameters
pca_perc = 90;
which_outcome = 'engel';
which_year = 1;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Initialize figure
figure
set(gcf,'position',[1 1 1400 450])
tiledlayout(1,3,"TileSpacing",'tight','padding','tight')

%% Histogram of outcomes

%% Outcome by laterality

%% Modeled probability of left by outcome for those who underwent surgery

% Remove patients without response
response = 'soz_lats';
empty_class = cellfun(@isempty,T.(response));
T(empty_class,:) = [];

% Do the model
just_spikes = 1; % ALl features
rm_non_temporal = 0; % All patients
combine_br = 1;

% Remove non temporal patients if desired
if rm_non_temporal == 1
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];
elseif rm_non_temporal == 2 % only include non-temporal (excludes diffuse and multifocal)
    extra = strcmp(T.soz_locs,'other cortex') | strcmp(T.soz_locs,'frontal');
    T(~extra,:) = [];
end

out =  classifier_wrapper(T,features,pca_perc,combine_br,just_spikes,rm_non_temporal);

assert(isequal(out.names,T.names))

% Find those who had left sided surgery
surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
left_surg = surg & strcmp(T.surg_lat,'left');
left_temporal_surg = surg & strcmp(T.surg_loc,'temporal');

% Get outcomes
outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
outcome = cellfun(@(x) parse_outcome_new(x,which_outcome),T.(outcome_name),'UniformOutput',false);
left_surg_good = left_surg & strcmp(outcome,'good');
left_surg_bad = left_surg & strcmp(outcome,'bad');

if 1
    
    oT = table(out.names,out.scores,out.class,out.all_pred,left_surg_good,left_surg_bad);
    oT(oT.Var5 == false & oT.Var6==false,:) = [];
    oT
end

% Plot
nexttile
unpaired_plot(out.scores(left_surg_good),out.scores(left_surg_bad),{'good','bad'},'Modeled probability of left')
title('Predicted probability for patients who underwent left sided surgery')

%% Model performance by good vs bad outcome???


end