function [T,rT] = laterality_features

%{
To do:zzzxz
- more error checking- really need to show basic components correct
- spatial normalization - figure this out. Also, decide how to deal with
negatives. V. concerning that model near identical w/ vs w/o normalization.
- try other classifier (random forest)
%}

which_atlas = 'aal';
which_outcome = 'ilae';
which_montage = 'car';
which_feature = 'ai';
restrict_mt = 1;
%which_thing = 'nelecs';
rm_no_surg = 1; % remove those who don't get surgery
randomize_lats = 0; % set to 1 to randomize soz laterality (null data)
do_sp_norm = 0; % spatially normalize

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

%% Initialize out file
if exist([plot_folder,'results.html'],'file') ~=0
    delete([plot_folder,'results.html'])
end
fid = fopen([plot_folder,'results.html'],'a');

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

%% get atlas
switch which_atlas
    case 'aal'
        atlas_names = data.aal_names;
        switch which_montage
            case 'car'
                atlas = data.all_aal;
            case 'bipolar'
                atlas = data.all_bipolar_aal;
        end
    case 'brainnetome'
        atlas_names = data.brainnetome_names;
        switch which_montage
            case 'car'
                atlas = data.all_brainnetome;   
            case 'bipolar'
                atlas = data.all_bipolar_brainnetome;
        end        
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
if rm_no_surg
    outcome(~resection_or_ablation) = nan; % make non resection or ablation nan
end

%% Info about outcome
fprintf(fid,'<b>Results</b>');
fprintf(fid,['<p>%d of %d patients had ablation or resection and available 2-year outcome data. Using 2-year'...
    ' %s outcome data, %d patients had a good outcome (defined as %s) and %d had poor outcome.'],...
    sum(resection_or_ablation&~isnan(outcome)),npts,which_outcome,sum(outcome==1),outcome_rule,sum(outcome==0));

%% Define laterality
old_bilat = strcmp(soz_lats,'bilateral') | strcmp(soz_lats,'diffuse');
unilat = strcmp(soz_lats,'left') | strcmp(soz_lats,'right');
bilat = nan(length(old_bilat),1);
bilat(old_bilat) = 1;
bilat(unilat) = 0;
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');

% Define laterality categories
nright = nansum(right_lat); nleft = nansum(left_lat); nbilat = nansum(bilat); nempty = sum(cellfun(@isempty,soz_lats));
assert(nright+nleft+nbilat+nempty == npts)
lat_num = nan(npts,1);
lat_num(cellfun(@isempty,soz_lats)) = 3;
lat_num(bilat==1) = 0;
lat_num(left_lat==1) = 1;
lat_num(right_lat==1) = 2;
assert(~any(isnan(lat_num)))

%% Randomize lats (null data)
% My coherence results do indeed turn non significant for the SOZ vs non
% SOZ laterality side
if randomize_lats
    new_lat_num = lat_num(randperm(npts));
    bilat = nan(length(old_bilat),1);
    bilat(new_lat_num == 0) = 1;
    bilat(ismember(new_lat_num,[1 2])) = 0;

    right_lat = zeros(npts,1);
    right_lat(new_lat_num==2) = 1;

    left_lat = zeros(npts,1);
    left_lat(new_lat_num==1) = 1;

    assert(nansum(bilat)==nbilat && sum(right_lat)==nright && sum(left_lat)==nleft)

end

%% Define L-R lateralizations
atlas_lat = lateralize_regions_simple(atlas_names);
elec_lats = cellfun(@(x) elec_broad(x,atlas_names,atlas_lat), atlas,'uniformoutput',false);

%% Get features
% Initialize table
T = table(outcome,bilat);
feat_names = {};

for which_thing = {'coh','bp','spikes','nelecs'}
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

    % Spatially normalize
    if sp_norm && do_sp_norm
        old_thing = thing;
        thing = spatial_normalization_fc(thing,locs);
        if 0
            tp = 110;
            figure
            nexttile
            turn_nans_gray(old_thing{tp}(:,:,6))
            nexttile
            turn_nans_gray(thing{tp}(:,:,6))
        end
    end

    %% Get feature
    switch which_feature
        case 'ai'
            % Get asymmetry index
            feat = get_ai(thing,last_dim,uni,elec_lats);

            if 0
                figure
                for ip = 50:npts
                    hold off
                    scatter3(locs{ip}(:,1),locs{ip}(:,2),locs{ip}(:,3),100,thing{ip},"filled")
                    hold on
                    text(locs{ip}(:,1),locs{ip}(:,2),locs{ip}(:,3),elec_lats{ip})
                    title(sprintf('%d AI: %1.2f',ip, feat(ip,1)))
                    pause
                end
            end
        case 'ai_90'
            feat = get_ai_prtile(thing,last_dim,uni,elec_lats);

            if 0
                figure
                for ip = 50:npts
                    hold off
                    scatter3(locs{ip}(:,1),locs{ip}(:,2),locs{ip}(:,3),100,thing{ip},"filled")
                    hold on
                    text(locs{ip}(:,1),locs{ip}(:,2),locs{ip}(:,3),elec_lats{ip})
                    title(sprintf('%d AI: %1.2f',ip, feat(ip,1)))
                    pause
                end
            end
        case 'sd'
            % weighted standard distance
            uw_feat = cellfun(@(x) weighted_standard_distance(x,[]),locs);
            
            if strcmp(which_thing,'nelecs')
                feat = uw_feat;
            else
                % first convert to nodal measure
                if uni == 0
                    nodal_thing = cellfun(@(x) squeeze(nanmean(x,2)),thing,'uniformoutput',false); % average across 2nd dimension
                else
                    nodal_thing = thing;
                end
                feat = nan(npts,last_dim); ...
                for ip = 1:npts
                    curr_locs = locs{ip};
                    curr_nodal_thing = nodal_thing{ip};
                    % take abs
                    curr_nodal_thing = abs(curr_nodal_thing);
                    for id = 1:last_dim
                        feat(ip,id) = weighted_standard_distance(curr_locs,curr_nodal_thing(:,id));
                    end
                end
                %feat = feat./uw_feat;
            end

            if 0
                figure
                for ip = 50:npts
                    hold off
                    scatter3(locs{ip}(:,1),locs{ip}(:,2),locs{ip}(:,3),100,abs(nodal_thing{ip}(:,1)),"filled")
                    hold on
                    text(locs{ip}(:,1),locs{ip}(:,2),locs{ip}(:,3),elec_lats{ip})
                    title(sprintf('%d SD: %1.2f',ip, feat(ip,1)))
                    pause
                end
            end
    end

    % prep table names
    tnames = cell(last_dim,1);
    for i = 1:last_dim
        tnames{i} = [which_thing{1},'_',num2str(i)];
    end
    feat_names = [feat_names;tnames];

    % Add features to table
    T = addvars(T,feat);
    T = splitvars(T,'feat','newVariableNames',tnames);
end

%% Compare asymmetry index between those with unilateral vs bilateral SOZ
if 0
    figure
    for f = 1:size(ai,2)
        nexttile
        unpaired_plot(feat(unilat==1,f),feat(bilat==1,f),{'unilateral','bilateral'},'thing');
    end
end

%% Compare asymmetry index between those with good vs bad outcome
if 0
    figure; set(gcf,'position',[100 100 300*size(feat,2) 350])
    tiledlayout(1,size(feat,2))
    for f = 1:size(feat,2)
        nexttile
        unpaired_plot(feat(outcome==1,f),feat(outcome==0,f),{'good outcome','bad outcome'},'thing');
    end
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
    print(gcf,[plot_folder,'feature_correlation'],'-dpng')
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
    print(gcf,[plot_folder,'univariate_bilateralty'],'-dpng')
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
    print(gcf,[plot_folder,'univariate_outcome'],'-dpng')
end

%% Try dimensionality reduction
% PCA
[coeff,score,latent,~,explained] = pca(all_feat,'rows','pairwise');
 %stem(cumsum(explained))

% first three explain almost 90% of the variance
ncoeffs = 3;
if 0
for f = 1:ncoeffs
    nexttile
    unpaired_plot(score(outcome==1,f),score(outcome==0,f),{'good','bad'},'thing');
end
end
coeff_names = cell(ncoeffs,1);
for i = 1:ncoeffs
    coeff_names{i} = ['coeff',num2str(i)];
end

rT = array2table([T.outcome,T.bilat,score(:,1:3)],'VariableNames',[{'bilat','outcome'},coeff_names']);

%% Do logistic regression
rng(0)
model = @(x,y,z) fitglm(x,'ResponseVar',y,'PredictorVars',z,'Distribution','binomial');
nb = 1e3;

% results in mean AUC of about 74% for full model
[~,~,AUC] = more_general_train_test(rT,nb,0.67,model,coeff_names,'outcome');

% Just elec num gets 65%  (null model)
[~,~,AUC_null] = more_general_train_test(T,nb,0.67,model,'nelecs_1','outcome');

% Just spikes gets 79%. So you really can't beat spikes. This is kind of
% cool. The presence of bilateral spikes really portends poor outcome
[trueClass,predClass,AUC_spikes] = more_general_train_test(T,nb,0.67,model,'spikes_1','outcome');



%% Confusion matrix
% Combine across nb
%{
tc = cell2mat(trueClass);
pc = cell2mat(predClass);
ac = nan(size(pc,1),1);
ac(pc>0.5) = 1; 
ac(pc<=0.5) = 0;
pc = ac;
any_nan = any(isnan([pc,tc]),2);
pc(any_nan) = [];
tc(any_nan) = [];
pc_cell = cell(size(pc,1),1);
tc_cell = cell(size(pc,1),1);
pc_cell(pc==1) = {'1'};
pc_cell(pc==0) = {'0'};
tc_cell(tc==1) = {'1'};
tc_cell(tc==0) = {'0'};

out = confusion_matrix(pc_cell,tc_cell,0);
%}

mean(AUC)
mean(AUC_null)
mean(AUC_spikes)

if 0
% What about predicting bilaterality
[~,~,AUC_bi] = more_general_train_test(rT,nb,0.67,model,coeff_names,'bilat');
[~,~,AUC_bi_null] = more_general_train_test(T,nb,0.67,model,'nelecs_1','bilat');
[~,~,AUC_bi_spikes] = more_general_train_test(T,nb,0.67,model,{'nelecs_1','spikes_1'},'bilat');

mean(AUC_bi)
mean(AUC_bi_null)
mean(AUC_bi_spikes)
end

end