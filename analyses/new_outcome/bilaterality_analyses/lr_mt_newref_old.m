function [T,features] =  lr_mt_newref
do_little_plots = 0;
do_big_plots = 0;
%which_outcome = 1; % engel = 1, ilae = 2
%which_outcome_year = 1;
which_sleep_stage = 3; % all = 1, wake =2, sleep = 3;


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
subplot_path = [plot_folder,'ai_subplots/'];
if ~exist(subplot_path,'dir')
    mkdir(subplot_path)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;

mt_data = load([inter_folder,'mt_out.mat']);
mt_data = mt_data.out;

%% get variables of interest
all_outcome = data.outcome; %outcome = all_outcome(:,which_outcome,which_outcome_year);
surgery = data.all_surgery;
soz_lats = data.all_soz_lats; 
soz_locs = data.all_soz_locs; 
names = data.all_names;
npts = length(names);
resection_lat = data.all_resec_lat;
ablation_lat = data.all_ablate_lat;
resection_loc = data.all_resec_loc;
ablation_loc = data.all_ablate_loc;
good_spikes = data.good_spikes;

%% All outcomes
engel_yr1 = all_outcome(:,1,1);
engel_yr2 = all_outcome(:,1,2);
ilae_yr1 = all_outcome(:,2,1);
ilae_yr2 = all_outcome(:,2,2);

%% Clean SOZ localizations and lateralities
soz_lats(strcmp(soz_lats,'diffuse')) = {'bilateral'}; % make diffuse be the same as bilateral
soz_locs(contains(soz_locs,'temporal')) = {'temporal'};

%% Consensus ablation or resection lat
surg_lat = cell(npts,1);
for i = 1:npts
    if strcmp(resection_lat{i},'na') && strcmp(ablation_lat{i},'na')
    elseif strcmp(resection_lat{i},'na')
        surg_lat{i} = ablation_lat{i};
    elseif strcmp(ablation_lat{i},'na')
        surg_lat{i} = resection_lat{i};
    else
        if ~strcmp(resection_lat{i},ablation_lat{i})
            error('what');
        end
        surg_lat{i} = resection_lat{i};
    end

end

%% Consensus ablation or reseciton loc
surg_loc = cell(npts,1);
for i = 1:npts
    if strcmp(resection_loc{i},'na') && strcmp(ablation_loc{i},'NA')
    elseif strcmp(resection_loc{i},'ATL')
        surg_loc{i} = 'temporal';
    elseif contains(ablation_loc{i},'temporal')
        surg_loc{i} = 'temporal';
    else
        surg_loc{i} = 'other';
    end
end


%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);
outcome(~resection_or_ablation) = {''}; % make non resection or ablation nan


%% Get features
% Initialize table
Ts = table(names,engel_yr1,engel_yr2,ilae_yr1,ilae_yr2,surgery,surg_lat,surg_loc,soz_locs,soz_lats);
feat_names_s = {};

for which_montage =[1 2 3] % machine,car, bipolar
    
    if which_montage == 1
        montage_text = 'machine';
    elseif which_montage == 2
        montage_text = 'car';
    elseif which_montage == 3
        montage_text = 'bipolar';
    end

    
    coh = mt_data.all_coh(:,which_montage,which_sleep_stage);
    pearson = mt_data.all_pearson(:,which_montage,which_sleep_stage);
    bp = mt_data.all_bp(:,which_montage,which_sleep_stage);
    rel_bp = mt_data.all_rel_bp(:,which_montage,which_sleep_stage);
    plv = mt_data.all_plv(:,which_montage,which_sleep_stage);
    rl = mt_data.all_rl(:,which_montage,which_sleep_stage);
    spikes = mt_data.all_spikes(:,which_montage,which_sleep_stage);
    xcor = mt_data.all_xcor(:,which_montage,which_sleep_stage);
    re = mt_data.all_re(:,which_montage,which_sleep_stage);
    se = mt_data.all_se(:,which_montage,which_sleep_stage);
    lags = mt_data.all_lags(:,which_montage,which_sleep_stage);

    % for now, make spike data nan if it was bad in the original
    % pipeline (I need to do my own validation eventually)
    %{
    for i = 1:npts
        %
        if good_spikes(i) == 0
            spikes{i} = nan(size(spikes{i}));
            rl{i} = nan(size(rl{i}));
        end
        %}
        %{
        if nanmean(spikes{i}) < 0.01
            spikes{i} = nan(size(spikes{i}));
            rl{i} = nan(size(rl{i}));
        end
        
    end
    %}
    
    % Loop over features
    for which_thing = {'spikes','rl','bp','se','pearson','xcor','coh','plv','re'}
        % Decide thing
        switch which_thing{1}
            case {'pearson','inter_pearson','near_pearson'}
                thing = pearson;
                uni = 0; % not univariate
                last_dim = 1; % one frequency
            case {'xcor','inter_xcor'}
                thing = xcor;
                uni = 0; % not univariate
                last_dim = 1; % one frequency
            case {'lags'}
                thing = lags;
                uni = 0; % not univariate
                last_dim = 1; % one frequency
            case {'coh','near_coh','inter_coh'}
                thing = coh;
                uni = 0; % not univariate
                last_dim = 6;%size(coh{1},3); % multiple frequencies
            case {'plv','near_plv','inter_plv'}
                thing = plv;
                uni = 0; % not univariate
                last_dim = 6;%size(coh{1},3); % multiple frequencies
            case {'re','inter_re'}
                thing = re;
                uni = 0; % not univariate
                last_dim = 6;%size(coh{1},3); % multiple frequencies
            case {'bp','inter_bp'}
                thing = bp;
                uni = 1; % univariate
                last_dim = 5;%size(bp{1},2); % multiple frequencies
            case {'rel_bp','inter_rel_bp'}
                thing = rel_bp;
                uni = 1; % univariate
                last_dim = 5;%size(bp{1},2); % multiple frequencies
            case {'se'}
                thing = se;
                uni = 1; % univariate
                last_dim = 1;%size(bp{1},2); % multiple frequencies
            case 'spikes'
                thing = spikes;
                uni = 1; % univariate
                last_dim = 1; % one frequency
            case 'nelecs'
                thing = cellfun(@(x) ones(length(x),1),spikes,'uniformoutput',false);
                last_dim = 1;
                uni = nan;
            case {'rl','inter_rl'}
                thing = rl;
                uni = 1;
                last_dim = 1;
        end
        
        labels = mt_data.all_labels(:,which_montage);
        

        %% Get asymmetry index
        %{
        %k = find(strcmp(names,'HUP167'));
        k = 100;
        calc_ai_sandbox(labels{k},thing{k},names{k},mt_data.all_labels{k,1},uni,last_dim,which_thing,subplot_path,do_little_plots)
        %}
        

        ai = cellfun(@(x,y,z,w) ...
            calc_ai_sandbox(x,y,z,w,uni,last_dim,which_thing,subplot_path,do_little_plots),...
            labels,thing,names,mt_data.all_labels(:,1),'uniformoutput',false);
    
        ai = cell2mat(ai);
    
    
        
        %% AI table
        tnames_s = cell(last_dim,1);
        for i = 1:last_dim
            tnames_s{i} = [which_thing{1},'_',num2str(i),'_',montage_text];
        end
        feat_names_s = [feat_names_s;tnames_s];
    
    
        Ts = addvars(Ts,ai);
        Ts = splitvars(Ts,'ai','newVariableNames',tnames_s);
    end

end

%% Remove redudnant features
%{
if sum(ismember(Ts.Properties.VariableNames,'spikes_1_bipolar'))>0
    Ts = removevars(Ts,'spikes_1_bipolar');
    feat_names_s(strcmp(feat_names_s,'spikes_1_bipolar')) = [];
end
if sum(ismember(Ts.Properties.VariableNames,'rl_1_bipolar'))>0
    Ts = removevars(Ts,'rl_1_bipolar');
    feat_names_s(strcmp(feat_names_s,'rl_1_bipolar')) = [];
end
%}

%% Pairwise correlations of all features
nfeatures = length(feat_names_s); % -2 to remove outcome and bilaterality
all_feat = table2array(Ts(:,size(Ts,2)-nfeatures+1:end));
feat_corr = corr(all_feat,'rows','pairwise','type','spearman');
if do_big_plots

    figure
    set(gcf,'position',[-300 78 1400 1200])
   % tiledlayout(4,4,'tilespacing','tight','Padding','tight')
    nexttile
    turn_nans_gray(feat_corr)
    xticks(1:nfeatures)
    xticklabels(cellfun(@(x) strrep(x,'_car',''),feat_names_s,'uniformoutput',false))
    yticks(1:nfeatures)
    yticklabels(cellfun(@(x) strrep(x,'_car',''),feat_names_s,'uniformoutput',false))
    caxis([-1 1])
    c = colorbar;
    c.Label.String = 'Correlation';
    %title('Correlation between L-R asymmetry indices')
    set(gca,'fontsize',15)

    figure
    set(gcf,'position',[-300 78 1400 1200])
   % tiledlayout(4,4,'tilespacing','tight','Padding','tight')

    cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

     
    for f = 1:length(feat_names_s)
        nexttile
        curr_feat = Ts.(feat_names_s{f});
        boxplot(Ts.(feat_names_s{f}),Ts.soz_lats,'colors',cols,'symbol','');
        hold on
        unique_lats = xticklabels;
        nlats = length(unique_lats);
        for il = 1:nlats
            curr_lats = strcmp(Ts.soz_lats,unique_lats{il});
            
            plot(il + randn(sum(curr_lats),1)*0.05,curr_feat(curr_lats),'o','color',cols(il,:))
        end
        %set(h,{'linew'},{2})
        %title(feat_names_s{f})
        ylabel(strrep(feat_names_s{f},'_car',''))
        yl = ylim;
        new_y = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];
        ylim(new_y)
        p = kruskalwallis(Ts.(feat_names_s{f}),Ts.soz_lats,'off');
        bon_p = 0.05/3;

        %% Only do post hoc tests if group comparison significant
        if p < 0.05
            % do post hoc
            lrp = ranksum(curr_feat(strcmp(Ts.soz_lats,'left')),curr_feat(strcmp(Ts.soz_lats,'right')));
            rbp = ranksum(curr_feat(strcmp(Ts.soz_lats,'right')),curr_feat(strcmp(Ts.soz_lats,'bilateral')));
            lbp = ranksum(curr_feat(strcmp(Ts.soz_lats,'left')),curr_feat(strcmp(Ts.soz_lats,'bilateral')));

            ybar1 = yl(1) + 1.06*(yl(2)-yl(1));
            ytext1 = yl(1) + 1.09*(yl(2)-yl(1));
            ybar2 = yl(1) + 1.18*(yl(2)-yl(1));
            ytext2 = yl(1) + 1.21*(yl(2)-yl(1));
            if lrp < bon_p
                plot([1 2],[ybar1 ybar1],'k-','linewidth',1)
                text(1.5,ytext1,get_asterisks_bonferroni(lrp,3),'horizontalalignment','center','fontsize',15)
            end
            if rbp < bon_p
                plot([2 3],[ybar1 ybar1],'k-','linewidth',1)
                text(2.5,ytext1,get_asterisks_bonferroni(rbp,3),'horizontalalignment','center','fontsize',15)
            end
            if lbp < bon_p
                plot([1 3],[ybar2 ybar2],'k-','linewidth',1)
                text(2,ytext2,get_asterisks_bonferroni(lbp,3),'horizontalalignment','center','fontsize',15)
            end
        end
        
        
       
        set(gca,'fontsize',15)
        
    end
    print(gcf,[plot_folder,'Fig2'],'-dpng')

end



%% Prep to output table, remove those with missing features
T = Ts(:,[1,size(Ts,2)-nfeatures+1:end]);
not_missing = ~any(ismissing(T),2);
T = Ts(not_missing,:);
features = feat_names_s;

%% Double check labels
if 0
good_labels = data.all_labels(not_missing,1);
for i = 1:length(good_labels)
    table(good_labels{i}(contains(good_labels{i},'A')|contains(good_labels{i},'B')|contains(good_labels{i},'C')))
    pause
end
end


end