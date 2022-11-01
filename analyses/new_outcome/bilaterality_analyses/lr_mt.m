function T =  lr_mt

which_outcome = 'ilae';
%which_montage = 'bipolar';

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
ilae = data.all_one_year_ilae;
engel = data.all_one_year_engel;
surgery = data.all_surgery;
soz_lats = data.all_soz_lats; soz_lats(strcmp(soz_lats,'diffuse')) = {'bilateral'};
npts = length(soz_lats);
names = data.all_names;
npts = length(names);
locs = data.all_native_locs;
bipolar_locs = data.all_native_bipolar_locs;



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
Ts = table(outcome,soz_lats);
feat_names_s = {};

for which_montage = {'car'}%;{'bipolar','car'}
    switch which_montage{1}
        case 'bipolar'
            coh = data.all_bipolar_coh;
            fc = data.all_bipolar_fc;
            bp = data.all_bp;
            spikes = data.all_spikes;
            labels = data.all_bipolar_labels;
        case 'car'  
            coh = data.all_coh;
            fc = data.all_fc;
            bp = data.all_bipolar_bp;
            spikes = data.all_spikes;
            labels = data.all_labels;
    end

    for which_thing = {'bp','spikes','coh','fc'}
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
            case 'near_coh'
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
        %[ai,signed] = cellfun(@(x,y) intra_mt_electrode_thing(x,y,uni,last_dim),labels,thing,'uniformoutput',false);
        %ai = cell2mat(ai);
        [signed,match] = cellfun(@(x,y) alt_intra_mt_electrode_thing(x,y,uni,last_dim,which_thing),labels,thing,'uniformoutput',false);
    
        signed = cell2mat(signed); ai = abs(signed);
        feat = ai;
    
        if 0
            figure
            tiledlayout(1,2)
            %unpaired_plot(signed(strcmp(soz_lats,'left')),signed(strcmp(soz_lats,'right')),{'left','right'},'signed diff')
            nexttile; boxplot(signed,soz_lats);
            nexttile; boxplot(ai,soz_lats);
        end
    
        
        %% Signed table
        tnames_s = cell(last_dim,1);
        for i = 1:last_dim
            tnames_s{i} = [which_thing{1},'_',num2str(i),'_',which_montage{1}];
        end
        feat_names_s = [feat_names_s;tnames_s];
    
    
        Ts = addvars(Ts,signed);
        Ts = splitvars(Ts,'signed','newVariableNames',tnames_s);
    end

end


if 0
    figure
    set(gcf,'position',[15 78 1400 350])
    tiledlayout(2,7,'tilespacing','tight','Padding','tight')
    for f = 1:size(Ts,2)-2
        nexttile
        %unpaired_plot(all_feat_s(strcmp(Ts.soz_lats,'left'),f),all_feat_s(strcmp(Ts.soz_lats,'right'),f),{'left','right'},feat_names_s{f});
        boxplot(Ts.(feat_names_s{f}),Ts.soz_lats)
        hold on
        ylabel(feat_names_s{f})
        p = kruskalwallis(Ts.(feat_names_s{f}),Ts.soz_lats,'off');
        yl = ylim;
        ybar = yl(1) + 1.05*(yl(2)-yl(1));
        ytext = yl(1) + 1.1*(yl(2)-yl(1));
        new_y = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
        plot([1 3],[ybar ybar],'k-','linewidth',2)
        text(2,ytext,sprintf('p = %1.3f',p),'horizontalalignment','center','fontsize',15)
        ylim(new_y)
        set(gca,'fontsize',15)
    end

end

%% Pairwise correlations of all features
nfeatures = size(Ts,2)-2; % -2 to remove outcome and bilaterality
all_feat = table2array(Ts(:,3:end));
feat_corr = corr(all_feat,'rows','pairwise');
if 0
    figure
    turn_nans_gray(feat_corr)
    xticks(1:nfeatures)
    xticklabels(feat_names_s)
    yticks(1:nfeatures)
    yticklabels(feat_names_s)
    colorbar
    title('Correlation between L-R asymmetry indices')
    set(gca,'fontsize',15)
    %print(gcf,[plot_folder,'feature_correlation'],'-dpng')
end


%no_nan  = ~isnan(Ts.fc_1 )& ~isnan(Ts.spikes_1) & ~isnan(Ts.bp_1);
%T = Ts(no_nan,:);
T = Ts(~any(ismissing(Ts),2),:);


end