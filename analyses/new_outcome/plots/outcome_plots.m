function outcome_plots

%% Parameters
which_year = 1;
outcome_approach = 'prob_comb';
which_model = 'spikes';

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

%% Initialize results file
fname = [plot_folder,'results.html'];
fid = fopen(fname,'a');
fprintf(fid,'<p><br><b>Concordance between spike-predicted laterality and surgical laterality is higher for patients with good surgical outcomes</b></br>');

fprintf(fid,['We examined the one-year surgical outcomes of patients who underwent resection '...
    'or laser ablation, studying both Engel and ILAE outcome classifications.']);
%% Load the model file
out = load([plot_folder,'ext_models.mat']);
out = out.all;

% which model
switch which_model
    case 'full'
        model = out.model(1).val(1);
    case 'spikes'
        model = out.model(2).val(1);
end

%% Run the mt_lr again just to get overall outcome stuff
T =  lr_mt(3);
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];
temporal_loc = contains(T.soz_locs,'temporal');
T(~temporal_loc,:) = [];
hup = contains(T.names,'HUP');
T(~hup,:) = [];

% check that outcomes match those in my manual validation table
if 0 
T(contains(T.surgery,'ablation')|contains(T.surgery,'Resection'),:)
% I double checked, they do. Note that some patients had "bilateral SOZs"
% but unilateral surgery. I count these in the outcomes and use the model
% concordant with the side of their surgery.
end


%% Initialize figure
figure
switch outcome_approach
    case 'prob_comb'
        set(gcf,'position',[1 1 1000 1000])
        tiledlayout(2,2,"TileSpacing",'compact','padding','tight')
    case 'auc_comb'
        set(gcf,'position',[1 1 1000 1000])
        tiledlayout(2,2,"TileSpacing",'compact','padding','tight')
    case 'prob_comb_alt'
        set(gcf,'position',[1 1 1000 1000])
        tiledlayout(2,2,"TileSpacing",'compact','padding','tight')
    otherwise
        set(gcf,'position',[1 1 1400 1000])
        tiledlayout(2,3,"TileSpacing",'compact','padding','tight')
end

% Prep stats for text
good_bad = nan(2,2); % engel, ilae; good, bad
prob_stats = nan(2,7); % engel, ilae; mean good, std good, mean bad, std bad, df, tstat, p

% Loop over outcome approaches
for io = 1:2
    if io == 1
        which_outcome = 'engel';
    elseif io == 2
        which_outcome = 'ilae';
    end


    % Anonymous function to define good outcome for first plot
    switch which_outcome
        case 'engel'
            good_outcome = @(x) strcmp(x(2),'A') | strcmp(x(2),'B') | strcmp(x(2),'C') | strcmp(x(2),'D');
            which_outcome_text = 'Engel';
            good_outcome_text = 'Engel I';
            bad_outcome_text = 'Engel II+';
        case 'ilae'
            good_outcome = @(x) contains(x,'2') | contains(x,'1');
            which_outcome_text = 'ILAE';
            good_outcome_text = 'ILAE 1-2';
            bad_outcome_text = 'ILAE 3+';
    
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


    
    %% Outcome analysis
    
    left = model.side(1).result;
    right = model.side(2).result;
    
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

    % how many good and bad outcome
    good_bad(io,1) = sum(good_outcome==1);
    good_bad(io,2) = sum(bad_outcome==1);
    
    
    switch outcome_approach
        case 'auc' % Get AUCs for good and bad outcomes
            % Hypothesis: if you had left surgery & good outcome, left model
            % AUC higher than if you had left surgery & bad outcome. Same for
            % right. So overall, the AUC for the corresponding model is higher
            % for good outcome.
            [XLg,YLg,~,AUCLg] = perfcurve(left.class(good_outcome),left.scores(good_outcome),...
                left.pos_class);
            [XLb,YLb,~,AUCLb] = perfcurve(left.class(bad_outcome),left.scores(bad_outcome),...
                left.pos_class);
            [XRg,YRg,~,AUCRg] = perfcurve(right.class(good_outcome),right.scores(good_outcome),...
                right.pos_class);
            [XRb,YRb,~,AUCRb] = perfcurve(right.class(bad_outcome),right.scores(bad_outcome),...
                right.pos_class);
    
            [pl,zl] = compare_independent_rocs(AUCLg,AUCLb,left.class(good_outcome),...
                left.class(bad_outcome),left.pos_class,left.pos_class);
            [pr,zr] = compare_independent_rocs(AUCRg,AUCRb,right.class(good_outcome),...
                right.class(bad_outcome),right.pos_class,right.pos_class);
    
            nexttile
            ll = plot(XLg,YLg,'linewidth',2);
            hold on
            lr = plot(XLb,YLb,':','linewidth',2);
            plot([0 1],[0 1],'k--','linewidth',2)
            xlabel('False positive rate')
            ylabel('True positive rate')
            legend([ll,lr],{sprintf('Left vs right/bilateral good: AUC = %1.2f',AUCLg),...
                sprintf('Left vs right/bilateral bad: AUC = %1.2f',AUCLb)},'fontsize',15,...
                'location','southeast')
            title(sprintf('Hanley-McNeil %s',get_p_text(pl)))
            set(gca,'fontsize',15)
    
            nexttile
            ll = plot(XRg,YRg,'linewidth',2);
            hold on
            lr = plot(XRb,YRb,':','linewidth',2);
            plot([0 1],[0 1],'k--','linewidth',2)
            xlabel('False positive rate')
            ylabel('True positive rate')
            legend([ll,lr],{sprintf('Right vs left/bilateral good: AUC = %1.2f',AUCRg),...
                sprintf('Right vs left/bilateral bad: AUC = %1.2f',AUCRb)},'fontsize',15,...
                'location','southeast')
            title(sprintf('Hanley-McNeil %s',get_p_text(pr)))
            set(gca,'fontsize',15)
        case 'auc_comb' % Get AUCs for good and bad outcomes
            % Hypothesis: if you had left surgery & good outcome, left model
            % AUC higher than if you had left surgery & bad outcome. Same for
            % right. So overall, the AUC for the corresponding model is higher
            % for good outcome.
            

    
        case 'prob' % hypothesis: the modeled probability of concordant laterality is higher for good outcome patients
            left_scores = left.scores;
            right_scores = right.scores;
    
            nexttile
            
            unpaired_plot(left_scores(good_outcome&left_surg),left_scores(bad_outcome&left_surg),...
                {'Good','Poor'},'Modeled probability of left laterality','para')
            %}
            %{
            unpaired_plot(left_scores(good_outcome&left_surg)-right_scores(good_outcome&left_surg),...
                left_scores(bad_outcome&left_surg)-right_scores(bad_outcome&left_surg),...
                {'Good','Poor'},'Modeled probability of left laterality','para')
            %}
            title({'Left-sided probability for patients','who underwent left-sided surgery'})
            xlim([0.5 2.5])
            set(gca,'fontsize',20)
    
            nexttile
            
            unpaired_plot(right_scores(good_outcome&right_surg),right_scores(bad_outcome&right_surg),{'Good','Poor'},...
                'Modeled probability of right laterality','para')
            %}
            %{
            unpaired_plot(right_scores(good_outcome&right_surg)-left_scores(good_outcome&right_surg),...
                right_scores(bad_outcome&right_surg)-left_scores(bad_outcome&right_surg),{'Good','Poor'},...
                'Modeled probability of right laterality','para')
            %}
            title({'Right-sided probability for patients','who underwent right-sided surgery'})
            xlim([0.5 2.5])
            set(gca,'fontsize',20)
    
        case 'prob_comb' % hypothesis: the modeled probability of concordant laterality is higher for good outcome patients
            left_scores = left.scores;
            right_scores = right.scores;
    
            nexttile
            % concordant lateralty: probability of left for those with left
            % surgery; probability of right for those with right surgery
            stats = unpaired_plot([left_scores(good_outcome&left_surg);right_scores(good_outcome&right_surg)],...
                [left_scores(bad_outcome&left_surg);right_scores(bad_outcome&right_surg)],...
                {good_outcome_text,bad_outcome_text},{'Modeled probability of','concordant laterality'},'para');
            set(gca().Children(3),'MarkerSize',10)
            set(gca().Children(4),'MarkerSize',10)
            title({'Surgery-model laterality concordance'})
            xlim([0.5 2.5])
            set(gca,'fontsize',20)

            % double check some stuff
            if 0
                table(T.names(left_surg),T.surg_lat(left_surg),T.ilae_yr1(left_surg),good_outcome(left_surg),bad_outcome(left_surg),left_scores(left_surg),right_scores(left_surg))
                table(T.names(right_surg),T.surg_lat(right_surg),T.ilae_yr1(right_surg),good_outcome(right_surg),bad_outcome(right_surg),left_scores(right_surg),right_scores(right_surg))
            end

            % engel, ilae; mean good, std good, mean bad, std bad, df, tstat, p
            prob_stats(io,:) = [stats.means(1) stats.sd(1) stats.means(2) stats.sd(2),...
                stats.df stats.tstat stats.p];

            % investigating the patients with high predicted concordant
            % laterality but poor outcome
            if 0
                table([T.names(bad_outcome&left_surg);T.names((bad_outcome&right_surg))],[left_scores(bad_outcome&left_surg);right_scores(bad_outcome&right_surg)])
            end
            % HUP138 (0.81) had a high modeled
            % probability of concordant laterality but poor outcome. They
            % had an ablation, and no repeat surgical evaluation. 

    
        case 'prob_comb_alt' % hypothesis: the modeled probability of concordant laterality is higher for good outcome patients
            % more complicated, takes difference in model scores. Slightly
            % more significant, but more complicated idea than above.
            left_scores = left.scores;
            right_scores = right.scores;
    
            nexttile
            % concordant lateralty: probability of left for those with left
            % surgery; probability of right for those with right surgery
            unpaired_plot([left_scores(good_outcome&left_surg)-right_scores(good_outcome&left_surg);...
                right_scores(good_outcome&right_surg)-left_scores(good_outcome&right_surg)],...
                [left_scores(bad_outcome&left_surg)-right_scores(bad_outcome&left_surg);...
                right_scores(bad_outcome&right_surg)-left_scores(bad_outcome&right_surg)],...
                {'Good','Poor'},{'Modeled probability of','concordant laterality'},'para')
            title({'Surgery-model laterality concordance'})
            xlim([0.5 2.5])
            set(gca,'fontsize',20)

            % double check some stuff
            if 0
                table(T.names(left_surg),T.surg_lat(left_surg),T.ilae_yr1(left_surg),good_outcome(left_surg),left_scores(left_surg),right_scores(left_surg))
                table(T.names(right_surg),T.surg_lat(right_surg),T.ilae_yr1(right_surg),good_outcome(right_surg),left_scores(right_surg),right_scores(right_surg))
            end
    
    end
end

%% Do text
fprintf(fid,[' %d of %d (%1.1f%%) had good one-year Engel outcomes (defined as Engel I), '...
    'and %d of %d (%1.1f%%) had poor Engel outcomes (defined as Engel 2+). '...
    '%d of %d (%1.1f%%) had good one-year ILAE outcomes (defined as ILAE 1-2), and '...
    '%d of %d (%1.1f%%) had poor ILAE outcomes (defined as ILAE 3+).'],good_bad(1,1),...
    sum(good_bad(1,:)),good_bad(1,1)/sum(good_bad(1,:))*100,...
    good_bad(1,2),...
    sum(good_bad(1,:)),good_bad(1,2)/sum(good_bad(1,:))*100,...
    good_bad(2,1),...
    sum(good_bad(2,:)),good_bad(2,1)/sum(good_bad(2,:))*100,...
    good_bad(2,2),...
    sum(good_bad(2,:)),good_bad(2,2)/sum(good_bad(2,:))*100);

fprintf(fid,[' We hypothesized that patients with a good surgical outcome would have a '...
    'higher modeled probability of SOZ laterality concordant with the side of surgery '...
    '(e.g., the modeled probability of having a left-sided SOZ would be higher in patients who '...
    ' had left-sided surgery and a good outcome than those who had left-sided surgery and a poor outcome). '...
    'We identified the model (incorporating spike rates as the only feature) corresponding to the side of surgery (e.g., for '...
    'patients who underwent left-sided surgery, we used the model predicting left vs. right/bilateral SOZ). '...
    'We compared the modeled probability of concordant laterality between patients with good and poor outcome. '...
    'There was a non-significant trend toward higher concordant model probabilities in '...
    'patients with good Engel outcomes (mean (SD) %1.2f (%1.2f)) than in patients with poor '...
    'Engel outcomes (%1.2f (%1.2f)) (<i>t</i>(%d) = %1.1f, %s). Mean concordant model probability '...
    'was significantly higher in patients with good ILAE outcomes (%1.2f (%1.2f)) than '...
    'in patients with poor ILAE outcomes (%1.2f (%1.2f)) (<i>t</i>(%d) = %1.1f, %s). '...
    'Together, these results suggest that a model trained to predict the SOZ using spike rate '...
    'asymmetry also predicts surgical outcome.</p>'],...
    prob_stats(1,1),prob_stats(1,2),prob_stats(1,3),prob_stats(1,4),prob_stats(1,5),prob_stats(1,6),...
    get_p_html(prob_stats(1,7)),...
    prob_stats(2,1),prob_stats(2,2),prob_stats(2,3),prob_stats(2,4),prob_stats(2,5),prob_stats(2,6),...
    get_p_html(prob_stats(2,7)));
% engel, ilae; mean good, std good, mean bad, std bad, df, tstat, p

%% Model performance by good vs bad outcome???
 print(gcf,[plot_folder,'Fig6'],'-dpng')


end

