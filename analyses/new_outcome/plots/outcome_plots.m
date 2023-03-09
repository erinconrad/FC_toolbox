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
tiledlayout(1,3,"TileSpacing",'tight','padding','compact')

%% Histogram of outcomes
% find those who had surgery
surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
outcome = T.(outcome_name); 
empty_outcome = cellfun(@isempty,outcome);
out_cat = categorical(outcome(surg&~empty_outcome));
cats = unique(out_cat);
nexttile
histogram(out_cat,cats)
hold on
yl = ylim;
yl_new = [yl(1) (yl(2)-yl(1))*1.3];
ybar = (yl(2)-yl(1))*1.1;
ytext = (yl(2)-yl(1))*1.2;
ylim(yl_new)
plot([1 3],[ybar ybar],'Color',[0.4660, 0.6740, 0.1880]	,'linewidth',2)
text(2,ytext,'Good outcome','fontsize',20,'HorizontalAlignment','center',...
    'color',[0.4660, 0.6740, 0.1880])
plot([4 8],[ybar ybar],'Color',[0.8500, 0.3250, 0.0980],'linewidth',2)
text(6,ytext,'Poor outcome','fontsize',20,'HorizontalAlignment','center',...
    'color',[0.8500, 0.3250, 0.0980])
plot([3.5,3.5],ylim, 'k--','linewidth',2)
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
    
    oT = table(out.names,out.scores,out.class,out.all_pred,left_surg_good,T.("spikes bipolar sleep"),...
        'VariableNames',{'name','prob','true','pred','good_out','Spike AI'});
    oT(left_surg_good == false & left_surg_bad==false,:) = [];
    oT
end

% Plot
nexttile
unpaired_plot(out.scores(left_surg_good),out.scores(left_surg_bad),{'Good','Poor'},'Modeled probability of left')
title({'Model-predicted concordance','for patients who underwent left sided surgery'})
set(gca,'FontSize',20)

%% Model performance by good vs bad outcome???


end