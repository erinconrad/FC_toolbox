function outcome_plots

%% Parameters
which_year = 1;
which_model = 'spikes';
which_refs = {'car','bipolar','machine'};


%% Get file locs
locations = fc_toolbox_locs;
plot_folder = locations.el_plots_folder;
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the file containing intermediate data
inter_folder = locations.el_data_folder;
mt_data = load([inter_folder,'mt_out_epilepsy_laterality.mat']);
mt_data = mt_data.out;

%% Load the MUSC outcome file
muscT = readtable([inter_folder,'LEN patient list research erin.xlsx']);


%% Initialize results file
for ir = 1:length(which_refs)

    if ir == 1
        fname = [plot_folder,'results.html'];
        fid = fopen(fname,'a');
        fprintf(fid,'<br><u><i>Concordance between spike-predicted laterality and surgical laterality is higher for patients with good surgical outcomes</i></u></br>');
        
        fprintf(fid,['We examined the one-year surgical outcomes of patients who underwent resection '...
            'or laser ablation, studying both Engel and ILAE outcome classifications.']);
    end
    
    %% Load the model file
    out = load([plot_folder,sprintf('ext_models_%s.mat',which_refs{ir})]);
    out = out.all;
    
    % which model
    switch which_model
        case 'full'
            model = out.approach(1).model(1).val(1);
        case 'spikes'
            model = out.approach(1).model(2).val(1);
    end
    
    %% Run the mt_lr again just to get overall outcome stuff
    T =  lr_mt(mt_data,3);
    empty_class = cellfun(@isempty,T.soz_lats);
    T(empty_class,:) = [];
    temporal_loc = contains(T.soz_locs,'temporal');
    T(~temporal_loc,:) = [];

    %% Figure out outcomes for HUP vs MUSC
    hup = contains(T.names,'HUP');
   % musc = contains(T.names,'MP');


    %T(~hup,:) = []; % Erin took this out 8/8

    %% Remove patients (should be one patient for bipolar) with nan feature
    features = T.Properties.VariableNames;
    spike_features = features(contains(features,'spikes') & contains(features,which_refs{ir}));
    nan_feature = isnan(T{:,spike_features});
    T(nan_feature,:) = [];    
    
    %% Initialize figure
    figure
    set(gcf,'position',[1 1 1400 1000])
        tiledlayout(2,3,"TileSpacing",'tight','padding','tight')
    
    % Prep stats for text
    good_bad = nan(2,2); % engel, ilae; good, bad
    prob_stats = nan(2,7); % engel, ilae; mean good, std good, mean bad, std bad, df, tstat, p
    auc_stats = nan(2,3); %engel, ilae; simple, complicated, delong p-value
    
    % Loop over outcome approaches (Engel vs ILAE), each one gets its own row
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
    
        %% A and D: Show overall outcomes
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
        text((1+sum(good))/2,ytext,'Good outcome','fontsize',15,'HorizontalAlignment','center',...
            'color',[0.4660, 0.6740, 0.1880])
        plot([sum(good)+1 length(good)],[ybar ybar],'Color',[0.8500, 0.3250, 0.0980],'linewidth',2)
        text((sum(good)+1+length(good))/2,ytext,'Poor outcome','fontsize',15,'HorizontalAlignment','center',...
            'color',[0.8500, 0.3250, 0.0980])
        plot([(sum(good)+sum(good)+1)/2,(sum(good)+sum(good)+1)/2],ylim, 'k--','linewidth',2)
        ylabel('Number of patients')
        title(sprintf('%s outcome',which_outcome_text))
        set(gca,'fontsize',15)
    
    
        
        %% B and E: See if modeled probability of concordant laterality is higher for good outcome patients
        % Get models
        left = model.side(1).result;
        right = model.side(2).result;
        
        % confirm that the patients all align with the outcome table
        assert(isequal(left.names,right.names))
        names = left.names;
        assert(isequal(T.names,names))
    
        % 224 did not have surgery, just planned to get it
        T.surgery(strcmp(T.names,'HUP224')) = {'none'};
        
        % Get some basic outcome stuff
        surg = (strcmp(T.surgery,'Laser ablation') | contains(T.surgery,'Resection'));
        outcome_bin = cellfun(@(x) parse_outcome_new(x,which_outcome),T.(outcome_name),'UniformOutput',false);
        good_outcome = strcmp(outcome_bin,'good') & surg == 1;
        bad_outcome = strcmp(outcome_bin,'bad') & surg == 1;
        left_surg = surg & strcmp(T.surg_lat,'left');
        right_surg = surg & strcmp(T.surg_lat,'right');
        npts = length(good_outcome);
    
        if 0
            table(T.names(surg),outcome_bin(surg),T.surgery(surg))
        end
    
        % Make sure no one had both left and right surg
        assert(sum(left_surg&right_surg)==0)
    
        % how many good and bad outcome
        good_bad(io,1) = sum(good_outcome==1);
        good_bad(io,2) = sum(bad_outcome==1);
        
        % Get the model scores
        left_scores = left.scores;
        right_scores = right.scores;
    
        % Get the concordant laterality scores for good and bad outcome
        concordant_lat_scores = nan(npts,1);
        concordant_lat_scores(left_surg == 1) = left_scores(left_surg==1);
        concordant_lat_scores(right_surg == 1) = right_scores(right_surg==1);
    
    
        nexttile
        % concordant lateralty: probability of left for those with left
        % surgery; probability of right for those with right surgery
        stats = unpaired_plot(concordant_lat_scores(good_outcome),concordant_lat_scores(bad_outcome),...
            {good_outcome_text,bad_outcome_text},{'Modeled probability of','concordant laterality'},'para');
        set(gca().Children(3),'MarkerSize',10)
        set(gca().Children(4),'MarkerSize',10)
        title({'Surgery-model laterality concordance'})
        xlim([0.5 2.5])
        set(gca,'fontsize',15)
    
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
        if 1
            table([T.names(bad_outcome&left_surg);T.names((bad_outcome&right_surg))],[left_scores(bad_outcome&left_surg);right_scores(bad_outcome&right_surg)])
        end
        % HUP138 (0.81) had a high modeled
        % probability of concordant laterality but poor outcome. They
        % had an ablation, and no repeat surgical evaluation. 
    
    
        %% C and F: Comparison between 5-Sense and 5-Sense + this model
        % Prep info about model
        scores = concordant_lat_scores;
        outcome_for_mdl = nan(length(scores),1);
        outcome_for_mdl(strcmp(outcome_bin,'good')) = 1; outcome_for_mdl(strcmp(outcome_bin,'bad')) = 0;
        nan_rows =  isnan(scores) | cellfun(@isempty,outcome_bin);
        non_nan_names = names(~nan_rows);
    
        % Get 5 sense score
        outT = five_sense_calculator(0,non_nan_names);
    
        % Reconcile 5-sense patients with my main table
        names_5sense = outT.names;
        prob_5sense = outT.prob;
        reconciled_prob = nan(length(names),1);
    
        for i = 1:length(names)
            % get row in the other outcome table corresponding to this name
            curr_name = names{i};
            row = strcmp(curr_name,names_5sense);
        
            if sum(row)~=1, continue; end
        
            reconciled_prob(i) = prob_5sense(row);
        end
    
        newT = table(names,outcome_bin,outcome_for_mdl,surg,left_surg,right_surg,reconciled_prob,left_scores,right_scores,scores);
        newT(isnan(newT.scores) | cellfun(@isempty,newT.outcome_bin),:) = []; % remove patients without model scores and patients without outcomes
        newT.outcome_for_mdl = logical(newT.outcome_for_mdl);
    
        % combined spike and 5-sense model to predict outcome
        probs = simple_loo_cv(newT,'outcome_for_mdl ~ reconciled_prob + scores','binomial');
        [X,Y,~,AUC] = perfcurve(newT.outcome_bin,probs,'good');
    
        % just 5-sense model to predict outcome
        simple_probs = simple_loo_cv(newT,'outcome_for_mdl ~ reconciled_prob','binomial');
        [simple_X,simple_Y,~,simple_AUC] = perfcurve(newT.outcome_bin,simple_probs,'good');
    
        % DeLong test comparing simple to complicated model
        pred = [probs,simple_probs];
        target = newT.outcome_for_mdl;
        % rm nans
        nan_rows = any(isnan(pred),2);
        pred(nan_rows,:) = [];
        target(nan_rows) = [];
    
        [ W ] = wilcoxon( pred, target );
        assert(abs(W(1)-AUC)<5e-3 && abs(W(2)-simple_AUC)<5e-3)
    
        [ S, S10, S01, V10, V01, theta ] = wilcoxonCovariance(pred,target);
        L = [1,-1]; % must always sum to 0
        alpha = 0.05; % significance level
        [ thetaP, thetaCI ] = wilcoxonConfidence(L, S, theta, alpha );
        % I should check that this gives the same answer as R -> I checked and
        % it did give the same answer as roc.test in R gives, and so I think
        % this is correct (although not significant)
    
        %writetable(table(target,pred),[plot_folder,'delong_test.csv'])
        auc_stats(io,:) = [simple_AUC,AUC,thetaP];
    
    
        %% do the plot
        nexttile
        cp = plot(X,Y,'linewidth',2);
        hold on
        sp = plot(simple_X,simple_Y,'linewidth',2);
        plot([0 1],[0 1],'k--','linewidth',2)
        xlabel('False positive rate')
        ylabel('True positive rate')
        legend([cp,sp],{sprintf('Model + pre-implant data: AUC = %1.2f',AUC),...
            sprintf('Pre-implant data only: AUC = %1.2f',simple_AUC)},...
            'location','southeast','fontsize',15)
        title('Outcome prediction')
        set(gca,'fontsize',15)
    
    
    end
    
    %% Do text
    if ir == 1
        fprintf(fid,[' %d of %d (%1.1f%%) patients had good one-year Engel outcomes (Engel I), '...
            'and %d of %d (%1.1f%%) had poor Engel outcomes (Engel 2+) (Fig. 4A). '...
            '%d of %d (%1.1f%%) patients had good one-year ILAE outcomes (ILAE 1-2), and '...
            '%d of %d (%1.1f%%) had poor ILAE outcomes (ILAE 3+) (Fig. 4C).'],good_bad(1,1),...
            sum(good_bad(1,:)),good_bad(1,1)/sum(good_bad(1,:))*100,...
            good_bad(1,2),...
            sum(good_bad(1,:)),good_bad(1,2)/sum(good_bad(1,:))*100,...
            good_bad(2,1),...
            sum(good_bad(2,:)),good_bad(2,1)/sum(good_bad(2,:))*100,...
            good_bad(2,2),...
            sum(good_bad(2,:)),good_bad(2,2)/sum(good_bad(2,:))*100);
        
        fprintf(fid,[' We hypothesized that patients with a good surgical outcome would have a '...
            'higher modeled probability of SOZ laterality concordant with the side of surgery. '...
            'We identified the spike rate model corresponding to the side of surgery (e.g., for '...
            'patients who underwent left-sided surgery, we used the model predicting left vs. right/bilateral SOZ). '...
            'Mean concordant model probability was significantly higher in patients with good Engel '...
            'outcomes (mean (SD) %1.2f (%1.2f)) than in patients with poor '...
            'Engel outcomes (%1.2f (%1.2f)) (<i>t</i>(%d) = %1.1f, %s) (Fig. 4B), and '...
            'win patients with good ILAE outcomes (%1.2f (%1.2f)) than '...
            'in patients with poor ILAE outcomes (%1.2f (%1.2f)) (<i>t</i>(%d) = %1.1f, %s) (Fig. 4E). '...
            'Together, these results suggest that a model trained to predict the SOZ using spike rate '...
            'asymmetry also predicts surgical outcome.</p>'],...
            prob_stats(1,1),prob_stats(1,2),prob_stats(1,3),prob_stats(1,4),prob_stats(1,5),prob_stats(1,6),...
            get_p_html(prob_stats(1,7)),...
            prob_stats(2,1),prob_stats(2,2),prob_stats(2,3),prob_stats(2,4),prob_stats(2,5),prob_stats(2,6),...
            get_p_html(prob_stats(2,7)));
        
        fprintf(fid,['<p>We also asked if spike rate asymmetry adds additional information to '...
            'predict surgical outcome beyond what is known using pre-implant data alone. '...
            'We used the recently published 5-SENSE score to estimate the probability of a single focal '...
            'seizure generator based on the pre-implant data. We used leave-one-patient-out cross-validation '...
            'to predict good vs. poor surgical outcome using a logistic regression model with the 5-SENSE score '...
            'as the only predictor. The 5-SENSE score predicted both Engel and ILAE outcomes (AUC = %1.2f '...
            'and %1.2f, respectively). We next applied the same cross-validation approach to a logistic '...
            'regression model incorporating both the 5-SENSE score and the modeled probability of concordant laterality. '...
            'The combined model had higher AUCs than the simple model for predicting both Engel and ILAE outcomes (AUC = '...
            '%1.2f and %1.2f, respectively), but neither combined model significantly outperformed '...
            'the 5-SENSE only model (DeLong test: %s for Engel outcome and %s for ILAE outcome) (Fig. 4C and F).'],...
            auc_stats(1,1),auc_stats(2,1),auc_stats(1,2),auc_stats(2,2),get_p_html(auc_stats(1,3)),get_p_html(auc_stats(2,3)));
        % engel, ilae; mean good, std good, mean bad, std bad, df, tstat, p

        fprintf(fid,[' Results were similar when we used spikes detected in bipolar and machine references to '...
            'build the SOZ laterality classifier (Fig. S4 and S5).</p>']);
    end
    
    %% Add subtitles
    annotation('textbox',[0 0.9 0.1 0.1],'String','A','LineStyle','none','fontsize',25)
    annotation('textbox',[0.33 0.9 0.1 0.1],'String','B','LineStyle','none','fontsize',25)
    annotation('textbox',[0.67 0.9 0.1 0.1],'String','C','LineStyle','none','fontsize',25)
    annotation('textbox',[0 0.4 0.1 0.1],'String','D','LineStyle','none','fontsize',25)
    annotation('textbox',[0.35 0.40 0.1 0.1],'String','E','LineStyle','none','fontsize',25)
    annotation('textbox',[0.67 0.4 0.1 0.1],'String','F','LineStyle','none','fontsize',25)
    
    
    
    if ir == 1
        print(gcf,[plot_folder,'Fig4'],'-dpng')
    elseif ir == 2
        print(gcf,[plot_folder,'FigS4'],'-dpng')
    elseif ir == 3
        print(gcf,[plot_folder,'FigS5'],'-dpng')
    end
    
    
end



end

