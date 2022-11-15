function SD = compare_features_laterality

%{

- How to do feature selection for binary variables?
- ALSO, THINK ABOUT BIPOLAR REFERENCE. Should be easy if just doing
left-right.
%}

%rng(0)

which_atlas = 'aal';
which_outcome = 'ilae';
restrict_mt = 0; % make anything outside of MT nans
spike_min = 0.1; % minimum spikes/min to calculate recruitment latency

coherence_freqs = {'delta','theta','alpha','beta','gamma','broadband'};

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

 %% get variables of interest
good_spikes = data.good_spikes;
%locs = data.all_locs;
locs = data.all_native_locs;
spike_rates = data.all_spikes;
labels = data.all_labels;
rl = data.all_rl;
fc = data.all_fc;
ilae = data.all_two_year_ilae;
engel = data.all_two_year_engel;
surgery = data.all_surgery;

anatomy = data.all_anatomy;
soz_lats = data.all_soz_lats;
bp = data.all_bp;
coh = data.all_coh;
bi_fc = data.all_bipolar_fc;
bi_coh = data.all_bipolar_coh;
bi_locs = data.all_bipolar_locs;


%% Get outcome
switch which_outcome
    case 'ilae'
        outcome = ilae;
    case 'engel'
        outcome = engel;
end

%% Find good and bad outcome
outcome_num = cellfun(@(x) parse_outcome(x,which_outcome),outcome);
outcome = outcome_num;

%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);

%% Make new variable for electrode nums
enums = cellfun(@(x) ones(length(x)),spike_rates,'uniformoutput',false);

%% get atlas
switch which_atlas
    case 'aal'
        atlas = data.all_aal;
        atlas_names = data.aal_names;
        bi_atlas = data.all_bipolar_aal;
    case 'brainnetome'
        atlas = data.all_brainnetome;
        atlas_names = data.brainnetome_names;
        bi_atlas = data.all_bipolar_brainnetome;
end

%% Define bilateral patients
old_bilat = strcmp(soz_lats,'bilateral') | strcmp(soz_lats,'diffuse');
unilat = strcmp(soz_lats,'left') | strcmp(soz_lats,'right');
bilat = nan(length(old_bilat),1);
bilat(old_bilat) = 1;
bilat(unilat) = 0;

%% Are all spikes nans if not good? - Yes
assert(isequal(cellfun(@(x) all(isnan(x)),spike_rates),~good_spikes))

%% Limit rl
% make all nans if spikes not good for that patient
rl = cellfun(@(x,y) rm_bad(x,y),rl,num2cell(good_spikes),'uniformoutput',false);

% Make those electrodes with few spikes nans
rl = cellfun(@(x,y) rm_few_spikes(x,y,spike_min),rl,spike_rates,'uniformoutput',false);

%% Define L-R lateralizations
if strcmp(which_atlas,'labels')
    elec_lats = cellfun(@(x,y) label_and_anatomy_lat_determination(x,y),labels,anatomy,'uniformoutput',false);
else
    atlas_lat = lateralize_regions_simple(atlas_names);
    elec_lats = cellfun(@(x) elec_broad(x,atlas_names,atlas_lat), atlas,'uniformoutput',false);
    elec_locs = cellfun(@(x) elec_broad(x,atlas_names,atlas_names), atlas,'uniformoutput',false);

    bi_elec_lats = cellfun(@(x) elec_broad(x,atlas_names,atlas_lat), bi_atlas,'uniformoutput',false);
    bi_elec_locs = cellfun(@(x) elec_broad(x,atlas_names,atlas_names), bi_atlas,'uniformoutput',false);
end

%% Restrict to mesial temporal?
% Get MT locs
mt = cellfun(@(x) find_mesial_temporal(x,which_atlas),elec_locs,'uniformoutput',false);

% Make anything that is not 1 in MT nans
if restrict_mt
    spike_rates = cellfun(@(x,y) assign_x_nan_where_y_0(x,y), spike_rates,mt,'uniformoutput',false);
    bp = cellfun(@(x,y) assign_x_nan_where_y_0(x,y), bp,mt,'uniformoutput',false);
    rl = cellfun(@(x,y) assign_x_nan_where_y_0(x,y), rl,mt,'uniformoutput',false);
    fc = cellfun(@(x,y) assign_x_nan_where_y_0(x,y), fc,mt,'uniformoutput',false);
    coh = cellfun(@(x,y) assign_x_nan_where_y_0(x,y), coh,mt,'uniformoutput',false);
    enums = cellfun(@(x,y) assign_x_nan_where_y_0(x,y), enums,mt,'uniformoutput',false);
    bi_fc = cellfun(@(x,y) assign_x_nan_where_y_0(x,y), bi_fc,mt,'uniformoutput',false);
    bi_coh = cellfun(@(x,y) assign_x_nan_where_y_0(x,y), bi_coh,mt,'uniformoutput',false);
end

%% Enum
mean_lr_num = cellfun(@(x,y) [nansum(x(strcmp(y,'L'))) nansum(x(strcmp(y,'R')))],...
    enums,elec_lats,'uniformoutput',false);
mean_lr_num = cell2mat(mean_lr_num);
enum_ai = asymmetry_index(mean_lr_num(:,1),mean_lr_num(:,2));
enuml = mean_lr_num(:,1);
enumr = mean_lr_num(:,2);

%% Get spike asymmetry index
% Get average left and right spike rates
mean_lr_spike_rate = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'))) nanmean(x(strcmp(y,'R')))],...
    spike_rates,elec_lats,'uniformoutput',false);
mean_lr_spike_rate = cell2mat(mean_lr_spike_rate);

% Asymmetry index
L = mean_lr_spike_rate(:,1);
R = mean_lr_spike_rate(:,2);
spike_ai = asymmetry_index(L,R);
spikesl = L;
spikesr = R;



%% RL asymmetry index
mean_lr_rl = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'))) nanmean(x(strcmp(y,'R')))],...
    rl,elec_lats,'uniformoutput',false);
mean_lr_rl = cell2mat(mean_lr_rl);
rl_ai = asymmetry_index(mean_lr_rl(:,1),mean_lr_rl(:,2));

%% BP asymmetry index
mean_lr_bp = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'),:),1)' nanmean(x(strcmp(y,'R'),:),1)'],...
    bp,elec_lats,'uniformoutput',false);
mean_lr_bp2 = cat(3,mean_lr_bp{:});
bp_ai = asymmetry_index(squeeze(mean_lr_bp2(:,1,:))',squeeze(mean_lr_bp2(:,2,:))');

%% FC
mean_lr_fc = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'),strcmp(y,'L')),'all') ...
nanmean(x(strcmp(y,'R'),strcmp(y,'R')),'all')],...
    fc,elec_lats,'uniformoutput',false);
mean_lr_fc = cell2mat(mean_lr_fc);
fc_ai = asymmetry_index(mean_lr_fc(:,1),mean_lr_fc(:,2));

%% Bi FC
mean_lr_bi_fc = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'),strcmp(y,'L')),'all') ...
nanmean(x(strcmp(y,'R'),strcmp(y,'R')),'all')],...
    bi_fc,bi_elec_lats,'uniformoutput',false);
mean_lr_bi_fc = cell2mat(mean_lr_bi_fc);
bi_fc_ai = asymmetry_index(mean_lr_bi_fc(:,1),mean_lr_bi_fc(:,2));

%% Coh
mean_lr_coh = cellfun(@(x,y) [squeeze(nanmean(x(strcmp(y,'L'),strcmp(y,'L'),:),[1 2])) squeeze(nanmean(x(strcmp(y,'R'),strcmp(y,'R'),:),[1 2]))],...
    coh,elec_lats,'uniformoutput',false);
mean_lr_coh2 = cat(3,mean_lr_coh{:});
coh_ai = asymmetry_index(squeeze(mean_lr_coh2(:,1,:))',squeeze(mean_lr_coh2(:,2,:))');

%% Bi Coh
mean_lr_bi_coh = cellfun(@(x,y) [squeeze(nanmean(x(strcmp(y,'L'),strcmp(y,'L'),:),[1 2])) squeeze(nanmean(x(strcmp(y,'R'),strcmp(y,'R'),:),[1 2]))],...
    bi_coh,bi_elec_lats,'uniformoutput',false);
mean_lr_bi_coh2 = cat(3,mean_lr_bi_coh{:});
bi_coh_ai = asymmetry_index(squeeze(mean_lr_bi_coh2(:,1,:))',squeeze(mean_lr_bi_coh2(:,2,:))');



%% Enum total
enum_tot = enuml+enumr;

%% Also get weighted dispersion of spikes
SD = cellfun(@(x,y) weighted_standard_distance(x,y),locs,spike_rates);
SD(SD>1e10) = nan;


% Also weighted dispersion of electrodes
SDE = cellfun(@(x) weighted_standard_distance(x,[]),locs);
SDE(SDE>1e10) = nan;

SDnorm = SD./SDE;

%T = table(bilat,enum_ai,spike_ai);
if 0
    figure; set(gcf,'Position',[100 100 1100 300])
    tiledlayout(1,4)
    nexttile
    unpaired_plot(fc_ai(bilat==0),fc_ai(bilat==1),{'Unilateral','Bilateral'},'Asymmetry index'); title('Pearson correlation-CAR')
    nexttile
    unpaired_plot(bi_fc_ai(bilat==0),bi_fc_ai(bilat==1),{'Unilateral','Bilateral'},'Asymmetry index'); title('Pearson correlation-bipolar')
    nexttile
    unpaired_plot(spike_ai(bilat==0),spike_ai(bilat==1),{'Unilateral','Bilateral'},'Asymmetry index'); title('Spikes')
    nexttile
    unpaired_plot(enum_ai(bilat==0),enum_ai(bilat==1),{'Unilateral','Bilateral'},'Asymmetry index'); title('Electrode number')
end

if 0
    figure; set(gcf,'Position',[100 100 1400 600])
    tiledlayout(2,6)
    for f = 1:6
        nexttile
        unpaired_plot(coh_ai(bilat==0,f),coh_ai(bilat==1,f),{'Unilateral','Bilateral'},'Asymmetry index');
        title(['CAR ',coherence_freqs{f}]);
    end
    for f = 1:6
        nexttile
        unpaired_plot(bi_coh_ai(bilat==0,f),bi_coh_ai(bilat==1,f),{'Unilateral','Bilateral'},'Asymmetry index');
        title(['Bipolar ',coherence_freqs{f}]);
    end
end

if 0
    figure; set(gcf,'Position',[100 100 1300 300])
    tiledlayout(1,5)
    for f = 1:5
        nexttile
        unpaired_plot(bp_ai(bilat==0,f),bp_ai(bilat==1,f),{'Unilateral','Bilateral'},'Asymmetry index');
        title(['CAR ',coherence_freqs{f}]);
    end

end




%% Outcome tables
if 0
    figure; set(gcf,'Position',[100 100 1100 300])
    tiledlayout(1,4)
    nexttile
    unpaired_plot(fc_ai(outcome==0),fc_ai(outcome==1),{'Bad','Good'},'Asymmetry index'); title('Pearson correlation-CAR')
    nexttile
    unpaired_plot(bi_fc_ai(outcome==0),bi_fc_ai(outcome==1),{'Bad','Good'},'Asymmetry index'); title('Pearson correlation-bipolar')
    nexttile
    unpaired_plot(spike_ai(outcome==0),spike_ai(outcome==1),{'Bad','Good'},'Asymmetry index'); title('Spikes')
    nexttile
    unpaired_plot(enum_ai(outcome==0),enum_ai(outcome==1),{'Bad','Good'},'Asymmetry index'); title('Electrode number')
end

if 0
figure
unpaired_plot(SD(outcome==0&resection_or_ablation),SD(outcome==1&resection_or_ablation),{'Bad','Good'},'Spatial dispersion'); title('Spikes')
end
%return

if 0
    figure; set(gcf,'Position',[100 100 1400 600])
    tiledlayout(2,6)
    for f = 1:6
        nexttile
        unpaired_plot(coh_ai(outcome==0,f),coh_ai(outcome==1,f),{'Bad','Good'},'Asymmetry index');
        title(['CAR ',coherence_freqs{f}]);
    end
    for f = 1:6
        nexttile
        unpaired_plot(bi_coh_ai(outcome==0,f),bi_coh_ai(outcome==1,f),{'Bad','Good'},'Asymmetry index');
        title(['Bipolar ',coherence_freqs{f}]);
    end
end

if 1
    figure; set(gcf,'Position',[100 100 1300 300])
    tiledlayout(1,5)
    for f = 1:5
        nexttile
        unpaired_plot(bp_ai(outcome==0,f),bp_ai(outcome==1,f),{'Bad','Good'},'Asymmetry index');
        title(['CAR ',coherence_freqs{f}]);
    end

end



%% Table
%T = table(outcome,bilat,enum_ai,spike_ai,fc_ai,bp_ai,coh_ai,SD,SDE,SDnorm,enuml,enumr,spikesl,spikesr);
%T = table(outcome,enuml,enumr,spikesl,spikesr);
T = table(outcome,enum_ai,bi_coh_ai,spike_ai);
T = splitvars(T);


%% Model to predict bilaterality
if  0
% Model parameters
N = 1e2;
perc_train = 2/3;
response = 'bilat';
model = @(x,y,z) fitglm(x,'ResponseVar',y,'PredictorVars',z,'Distribution','binomial');

% Define predictors
null_predictors = {'spike_ai'};
full_predictors = {'bp_ai_1'};

% Remove rows with missing predictors or response
missing_response = ismissing(T.(response));
missing_predictors = cellfun(@(x) ismissing(T.(x)),full_predictors,'UniformOutput',false);
missing_predictors = horzcat(missing_predictors{:});
missing_anything = any([missing_predictors,missing_response],2);
T(missing_anything,:) = [];

% Null model
[trueClass_null,predClass_null,AUC_null] = more_general_train_test(T,N,perc_train,model,null_predictors,response);


% FUll model
[trueClass,predClass,AUC] = more_general_train_test(T,N,perc_train,model,full_predictors,response);
end



%% Model to predict outcome
% Model parameters
N = 1e2;
perc_train = 2/3;
response = 'outcome';
model = @(x,y,z) fitglm(x,'ResponseVar',y,'PredictorVars',z,'Distribution','binomial');

% Define predictors
null_predictors = {'enum_ai'};
full_predictors = {'enum_ai','bi_coh_ai_6'};

% Remove rows with missing predictors or response
missing_response = ismissing(T.(response));
missing_predictors = cellfun(@(x) ismissing(T.(x)),full_predictors,'UniformOutput',false);
missing_predictors = horzcat(missing_predictors{:});
missing_anything = any([missing_predictors,missing_response],2);

% also remove those without resection/ablation
T(missing_anything | ~resection_or_ablation,:) = [];
%T(~resection_or_ablation,:) = [];

% Null model
[trueClass_null,predClass_null,AUC_null] = more_general_train_test(T,N,perc_train,model,null_predictors,response);


% FUll model
[trueClass,predClass,AUC] = more_general_train_test(T,N,perc_train,model,full_predictors,response);

nanmean(AUC_null)
nanmean(AUC)


end

function x = assign_x_nan_where_y_0(x,y)

if ndims(x) == 3
    x(y==0,:,:) = nan;
    x(:,y==0,:) = nan;

elseif size(x,2) > 1
    if size(x,2) == size(x,1)
        x(y==0,:) = nan;
        x(:,y==0) = nan;
    else
        x(y==0,:) = nan;
    end
else
    x(y==0) = nan;
end

end

function x = rm_bad(x,y)

if y == 0
x(:) = nan;
end

end

function x = rm_few_spikes(x,y,min_rate)

x(y<min_rate) = nan;

end

