function outcome_plots(T,features)

%% Parameters
pca_perc = 95;
which_outcome = 'engel';
which_year = 1;
outcome_approach = 'prob';
only_temp = 0;
direct_model = 0;

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
tiledlayout(1,3,"TileSpacing",'tight','padding','compact')

%% Define good outcome
switch which_outcome
    case 'engel'
        good_outcome = @(x) strcmp(x(2),'A') | strcmp(x(2),'B') | strcmp(x(2),'C') | strcmp(x(2),'D');
    case 'ilae'
        good_outcome = @(x) contains(x,'2') | contains(x,'1');

end

%% Histogram of outcomes
% find those who had surgery
surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
outcome = T.(outcome_name); 
empty_outcome = cellfun(@isempty,outcome);
out_cat = categorical(outcome(surg&~empty_outcome));
cats = unique(out_cat);
good = arrayfun(@(x) good_outcome(char(x)),cats);

nexttile
histogram(out_cat,cats)
hold on
yl = ylim;
yl_new = [yl(1) (yl(2)-yl(1))*1.3];
ybar = (yl(2)-yl(1))*1.1;
ytext = (yl(2)-yl(1))*1.2;
ylim(yl_new)
plot([1 sum(good)],[ybar ybar],'Color',[0.4660, 0.6740, 0.1880]	,'linewidth',2)
text((1+sum(good))/2,ytext,'Good outcome','fontsize',20,'HorizontalAlignment','center',...
    'color',[0.4660, 0.6740, 0.1880])
plot([sum(good)+1 length(good)],[ybar ybar],'Color',[0.8500, 0.3250, 0.0980],'linewidth',2)
text((sum(good)+1+length(good))/2,ytext,'Poor outcome','fontsize',20,'HorizontalAlignment','center',...
    'color',[0.8500, 0.3250, 0.0980])
plot([(sum(good)+sum(good)+1)/2,(sum(good)+sum(good)+1)/2],ylim, 'k--','linewidth',2)
ylabel('Number of patients')
title('Engel outcome')
set(gca,'fontsize',20)


%% Outcome by laterality
left_surg = surg & strcmp(T.surg_lat,'left');
right_surg = surg & strcmp(T.surg_lat,'right');
out_cat_left = categorical(outcome(left_surg&~empty_outcome));
out_cat_right = categorical(outcome(right_surg&~empty_outcome));
nexttile
histogram(out_cat_left,cats,'facecolor','none','edgecolor',[0, 0.4470, 0.7410],...
    'linewidth',2)
hold on
histogram(out_cat_right,cats,'facecolor','none','edgecolor',[0.8500, 0.3250, 0.0980],...
    'linewidth',2)
legend({'Left','Right'},'location','northeast','fontsize',20)
ylabel('Number of patients')
title('Engel outcome by surgery laterality')
set(gca,'fontsize',20)


%% Modeled probability of left by outcome for those who underwent surgery
if direct_model



    
    % Just predict outcome
    % Get outcomes
    outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
    outcome = cellfun(@(x) parse_outcome_new(x,which_outcome),T.(outcome_name),'UniformOutput',false);

    % Find those who had left sided surgery
    surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));

    % Only want outcome for those who had surgery
    outcome(~surg) = {''};
    T.outcome = outcome;

    % proof of concept
    if 1
        
        good = strcmp(outcome,'good');
        left = strcmp(T.soz_lats,'left');
        right = strcmp(T.soz_lats,'right');
        curr = T.('spikes machine sleep');
        figure
        nexttile
        unpaired_plot(abs(curr(good)),abs(curr(~good)),{'good','bad'},'spikes')
        nexttile
        unpaired_plot(curr(left),curr(right),{'left','right'},'spikes')
    end

    % remove empty class
    empty_class = cellfun(@isempty,T.outcome);
    T(empty_class,:) = [];

    % model
    just_spikes = 1; % ALl features
    rm_non_temporal = 0; % All patients
    combine_br = 0;
    out =  classifier_wrapper(T,features,pca_perc,combine_br,just_spikes,rm_non_temporal,'outcome');
    [X,Y,~,AUC] = perfcurve(out.class,out.scores,out.pos_class);


    nexttile
    lg = plot(X,Y,'linewidth',2);
    hold on
    plot([0 1],[0 1],'k--','linewidth',2)
    xlabel('False positive rate')
    ylabel('True positive rate')
    legend(sprintf('AUC = %1.2f',AUC),'fontsize',20,...
        'location','southeast')
    title({'Model performance by outcome'})
    set(gca,'fontsize',20)


else

%{
Approaches to this analysis:
- Modeled probability of left by good vs bad outcome (no difference)
- Outcome (continuous variable) according to accurate vs inaccurate
- AUC according to good and bad outcome
%} 

% Remove patients without response
response = 'soz_lats';
empty_class = cellfun(@isempty,T.(response));
T(empty_class,:) = [];

% Do the model
just_spikes = 1; % ALl features
rm_non_temporal = 0; % All patients

% Remove non temporal patients if desired
if rm_non_temporal == 1
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];
elseif rm_non_temporal == 2 % only include non-temporal (excludes diffuse and multifocal)
    extra = strcmp(T.soz_locs,'other cortex') | strcmp(T.soz_locs,'frontal');
    T(~extra,:) = [];
end

% Do the models
left =  classifier_wrapper(T,features,pca_perc,1,just_spikes,rm_non_temporal,[]);
right =  classifier_wrapper(T,features,pca_perc,2,just_spikes,rm_non_temporal,[]);

% Ensure that patient order in table lines up with that of model output
assert(all(strcmp(T.names,left.names))); assert(all(strcmp(T.names,right.names))) 

% Find those who had surgery
if only_temp
    surg = ((strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'))) & strcmp(T.surg_loc,'temporal');
else
    surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
end
left_surg = surg & strcmp(T.surg_lat,'left');
right_surg = surg & strcmp(T.surg_lat,'right');

% Match the score to the surg
scores = nan(length(left_surg),1);
left_scores = left.scores;
right_scores = right.scores;
scores(left_surg) = left_scores(left_surg);
scores(right_surg) = right_scores(right_surg);

% Get outcomes
outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
outcome = cellfun(@(x) parse_outcome_new(x,which_outcome),T.(outcome_name),'UniformOutput',false);
surg_good = surg & strcmp(outcome,'good');
surg_bad = surg & strcmp(outcome,'bad');

%outcome_num =cellfun(@(x)  parse_outcome_num(x,which_outcome), T.(outcome_name));

if 0
    
    table(T.names(surg_good|surg_bad),scores(surg_good|surg_bad),surg_good(surg_good|surg_bad))
    table(T.names((surg_good|surg_bad)&left_surg),scores((surg_good|surg_bad)&left_surg),surg_good((surg_good|surg_bad)&left_surg))
end

% Plot
nexttile

switch outcome_approach
    case 'prob'
        unpaired_plot(scores(surg_good),scores(surg_bad),{'Good','Poor'},'Modeled probability of concordant laterality')
        %unpaired_plot(left_scores(surg_good&left_surg),left_scores(surg_bad&left_surg),{'Good','Poor'},'Modeled probability of concordant laterality')
        title({'Prediction concordance for patients','who underwent surgery'})
        xlim([0.5 2.5])
        set(gca,'fontsize',20)
    case 'auc'
        
        [Xg,Yg,~,AUCg] = perfcurve(left.class(surg_good&left_surg),left.scores(surg_good&left_surg),left.pos_class);
        [Xb,Yb,~,AUCb] = perfcurve(left.class(surg_bad&left_surg),left.scores(surg_bad&left_surg),left.pos_class);

        lg = plot(Xg,Yg,'linewidth',2);
        hold on
        lb = plot(Xb,Yb,':','linewidth',2);
        plot([0 1],[0 1],'k--','linewidth',2)
        xlabel('False positive rate')
        ylabel('True positive rate')
        legend([lg,lb],{sprintf('Good outcome: AUC = %1.2f',AUCg),...
            sprintf('Bad outcome: AUC = %1.2f',AUCb)},'fontsize',20,...
            'location','southeast')
        title({'Model performance by outcome'})
        set(gca,'fontsize',20)

        
end

%% Model performance by good vs bad outcome???
print(gcf,[plot_folder,'Fig6'],'-dpng')


end

end