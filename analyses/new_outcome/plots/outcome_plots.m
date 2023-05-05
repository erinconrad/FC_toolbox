function outcome_plots

%% Parameters
which_outcome = 'engel';
which_year = 2;
outcome_approach = 'prob';
which_model = 'full';

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

%% Load the model file
out = load([plot_folder,'models.mat']);
out = out.out;

%% Run the mt_lr again just to get overall outcome stuff
[T,features] =  lr_mt(3);
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];
temporal_loc = contains(T.soz_locs,'temporal');
T(~temporal_loc,:) = [];

%% Initialize figure
figure
set(gcf,'position',[1 1 1400 450])
tiledlayout(1,3,"TileSpacing",'tight','padding','loose')

%% Anonymous function to define good outcome for first two plots
switch which_outcome
    case 'engel'
        good_outcome = @(x) strcmp(x(2),'A') | strcmp(x(2),'B') | strcmp(x(2),'C') | strcmp(x(2),'D');
        which_outcome_text = 'Engel';
    case 'ilae'
        good_outcome = @(x) contains(x,'2') | contains(x,'1');
        which_outcome_text = 'ILAE';

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
title(sprintf('%s outcome',which_outcome_text))
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

%% Outcome analysis
nexttile
% Get some basic model stuff
switch which_model
    case 'full'
        left = out.full_model.left;
        right = out.full_model.right;
    case 'spikes'
        left = out.spike_model.left;
        right = out.spike_model.right;

end

assert(isequal(left.names,right.names))
names = left.names;
assert(isequal(T.names,names))

% Get some basic outcome stuff
surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
outcome_bin = cellfun(@(x) parse_outcome_new(x,which_outcome),T.(outcome_name),'UniformOutput',false);
good_outcome = strcmp(outcome_bin,'good') & surg == 1;
bad_outcome = strcmp(outcome_bin,'bad') & surg == 1;
left_surg = surg & strcmp(T.surg_lat,'left');
right_surg = surg & strcmp(T.surg_lat,'right');
npts = length(good_outcome);


switch outcome_approach
    case 'auc' % Get AUCs for good and bad outcomes
        % Hypothesis: if you had left surgery & good outcome, left model
        % AUC higher than if you had left surgery & bad outcome. Same for
        % right. So overall, the AUC for the corresponding model is higher
        % for good outcome.
        [~,~,~,AUCLg] = perfcurve(left.class(good_outcome),left.scores(good_outcome),...
            left.pos_class);
        [~,~,~,AUCLb] = perfcurve(left.class(bad_outcome),left.scores(bad_outcome),...
            left.pos_class);
        [~,~,~,AUCRg] = perfcurve(right.class(good_outcome),right.scores(good_outcome),...
            right.pos_class);
        [~,~,~,AUCRb] = perfcurve(right.class(bad_outcome),right.scores(bad_outcome),...
            right.pos_class);
    case 'prob' % hypothesis: the modeled probability of concordant laterality is higher for good outcome patients
        scores = nan(npts,1);
        left_scores = left.scores;
        right_scores = right.scores;
        scores(left_surg) = left_scores(left_surg);
        scores(right_surg) = right_scores(right_surg);

        unpaired_plot(scores(good_outcome),scores(bad_outcome),{'Good','Poor'},'Modeled probability of concordant laterality')
        title({'Prediction concordance for patients','who underwent surgery'})
        xlim([0.5 2.5])
        set(gca,'fontsize',20)

end


%% Model performance by good vs bad outcome???
print(gcf,[plot_folder,'Fig6'],'-dpng')


end

