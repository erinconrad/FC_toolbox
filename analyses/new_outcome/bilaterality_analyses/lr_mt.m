function lr_mt

which_outcome = 'ilae';
which_montage = 'car';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;


%% get variables of interest
ilae = data.all_two_year_ilae;
engel = data.all_two_year_engel;
surgery = data.all_surgery;
soz_lats = data.all_soz_lats;
npts = length(soz_lats);
names = data.all_names;
npts = length(names);
locs = data.all_native_locs;
bipolar_locs = data.all_native_bipolar_locs;

switch which_montage
    case 'bipolar'
        coh = data.all_bipolar_coh;
        fc = data.all_bipolar_fc;
        bp = data.all_bp;
        spikes = data.all_spikes;
        locs = bipolar_locs;
        labels = data.all_bipolar_labels;
    case 'car'  
        coh = data.all_coh;
        fc = data.all_fc;
        bp = data.all_bipolar_bp;
        spikes = data.all_spikes;
        labels = data.all_labels;
end

%% Get outcome
switch which_outcome
    case 'ilae'
        outcome = ilae;
    case 'engel'
        outcome = engel;
end

%% Find good and bad outcome
outcome_num = cellfun(@(x) parse_outcome(x,which_outcome),outcome);
[~,outcome_rule] = parse_outcome('',which_outcome);
outcome = outcome_num;

%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);
outcome(~resection_or_ablation) = nan; % make non resection or ablation nan

%% Define laterality
old_bilat = strcmp(soz_lats,'bilateral') | strcmp(soz_lats,'diffuse');
unilat = strcmp(soz_lats,'left') | strcmp(soz_lats,'right');
bilat = nan(length(old_bilat),1);
bilat(old_bilat) = 1;
bilat(unilat) = 0;
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');

%% Get features
% Initialize table
T = table(outcome,bilat);
feat_names = {};

Ts = table(soz_lats);
feat_names_s = {};

for which_thing = {'spikes','coh','bp'}
    % Decide thing
    switch which_thing{1}
        case 'fc'
            thing = fc;
            uni = 0;
            last_dim = 1;
            sp_norm = 1;
        case 'coh'
            thing = coh;
            uni = 0;
            last_dim = size(coh{1},3);
            sp_norm = 1;
        case 'bp'
            thing = bp;
            uni = 1;
            last_dim = size(bp{1},2);
            sp_norm = 0;
        case 'spikes'
            thing = spikes;
            uni = 1;
            last_dim = 1;
            sp_norm = 0;

        case 'nelecs'
            thing = cellfun(@(x) ones(length(x),1),spikes,'uniformoutput',false);
            last_dim = 1;
            uni = nan;
            sp_norm = 0;
    end

    %% Get intra
    [ai,signed] = cellfun(@(x,y) intra_mt_electrode_thing(x,y,uni,last_dim),labels,thing,'uniformoutput',false);
    ai = cell2mat(ai);
    signed = cell2mat(signed);
    feat = ai;

    %% Main table
     % prep table names
    tnames = cell(last_dim,1);
    for i = 1:last_dim
        tnames{i} = [which_thing{1},'_',num2str(i)];
    end
    feat_names = [feat_names;tnames];

    % Add features to table
    T = addvars(T,feat);
    T = splitvars(T,'feat','newVariableNames',tnames);

    %% Signed table
    tnames_s = cell(last_dim,1);
    for i = 1:last_dim
        tnames_s{i} = [which_thing{1},'_',num2str(i)];
    end
    feat_names_s = [feat_names_s;tnames_s];


    Ts = addvars(Ts,signed);
    Ts = splitvars(Ts,'signed','newVariableNames',tnames_s);

end

%% Pairwise correlations of all features
nfeatures = size(T,2)-2; % -2 to remove outcome and bilaterality
all_feat = table2array(T(:,3:end));
feat_corr = corr(all_feat,'rows','pairwise');
if 0
    figure
    turn_nans_gray(feat_corr)
    xticks(1:nfeatures)
    xticklabels(feat_names)
    yticks(1:nfeatures)
    yticklabels(feat_names)
    colorbar
    title('Correlation between L-R asymmetry indices')
    set(gca,'fontsize',15)
    %print(gcf,[plot_folder,'feature_correlation'],'-dpng')
end

%% Univariate analysis of each feature with bilaterality
if 0
    figure
    set(gcf,'position',[15 78 1377 719])
    tiledlayout(3,5,'tilespacing','tight','Padding','tight')
    for f = 1:nfeatures
        nexttile
        unpaired_plot(all_feat(T.bilat==0,f),all_feat(T.bilat==1,f),{'unilateral','bilateral'},feat_names{f});
    end
    %print(gcf,[plot_folder,'univariate_bilateralty'],'-dpng')
end

%% Univariate analysis of each feature with laterality
if 0
    figure
    set(gcf,'position',[15 78 1377 719])
    tiledlayout(3,5,'tilespacing','tight','Padding','tight')
    for f = 1:size(Ts,2)-1
        nexttile
        unpaired_plot(all_feat(strcmp(Ts.soz_lats,'left'),f),all_feat(strcmp(Ts.soz_lats,'right'),f),{'left','right'},feat_names{f});

    end
end

%% Univariate analysis of each feature with outcome
if 1
    figure
    set(gcf,'position',[15 78 1377 719])
    tiledlayout(3,5,'tilespacing','tight','Padding','tight')
    for f = 1:nfeatures
        nexttile
        unpaired_plot(all_feat(T.outcome==1,f),all_feat(T.outcome==0,f),{'good','bad'},feat_names{f});
    end
   % print(gcf,[plot_folder,'univariate_outcome'],'-dpng')
end


end