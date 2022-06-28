function T = outcome_table(outcome,predictors)

%{
Inputs: 
- outcome: 1 or 0 indicating good or bad outcome
- predictor struct: an npredictor x 1 structure array containing the
following fields
   - type: a string indicating the type of predictor. Acceptable
   options are "binary","categorical",or "continuous"
    - X: an Nx1 array or cell array, where N is the same length as
   outcome
   - name: a string indicting the name of the predictor


% Need to figure out a 2x3 fisher table for categorical
%}

%% Get number of predictors
np = length(predictors);

%% Prep cell array representing table, start empty
outcome_cell = {};

%% Get the identity of nan outcomes
nan_outcome = isnan(outcome);
npts = sum(~nan_outcome);

%% Loop over predictors
for ip = 1:np
    
 
    % Get type and thing
    type = predictors(ip).type;
    X = predictors(ip).X;
    Xname = predictors(ip).name;
    
    
    %% if the type is continuous, just single row
    if strcmp(type,'continuous')
        
        % Get variable for two outcome groups
        X1 = X(outcome==1);
        X0 = X(outcome==0);
        
        % median and IQR
        medianX1 = nanmedian(X1);
        iqrX1 = prctile(X1,[25 75]);
        medianX0 = nanmedian(X0);
        iqrX0 = prctile(X0,[25 75]);
        
        % Do rank sum
        [p,~,stats] = ranksum(X1,X0);
        Tpos = stats.ranksum;
        
        %% Make row in cell
        outcome_cell = [outcome_cell;{Xname},...
            {sprintf('%1.1f (%1.1f-%1.1f)',medianX1,iqrX1(1),iqrX1(2))},...
            {sprintf('%1.1f (%1.1f-%1.1f)',medianX0,iqrX0(1),iqrX0(2))},...
            {Tpos},...
            {sprintf('%s',simple_p_text(p))}];
            
        
        
    end
    
    %% Categorical
    if strcmp(type,'categorical')
        error('you cannot do this yet')
    end
    
    %% Binary
    if strcmp(type,'binary')
        
        % Do cross tab to get 2x2 table and labels
        [tab2,~,~,labels] = crosstab(outcome,X);
        cat_labels = labels(:,2);
        prct_tab2 = tab2/npts;
        
        % do fisher exact test
        [~,p,stats] = fishertest(tab2);
        or = stats.OddsRatio;
        ci = stats.ConfidenceInterval;
        
        % Numbers in each category
        ncats = length(cat_labels);
    
        % First row just has name
        outcome_cell = [outcome_cell;{Xname},{''},{''},{''},{''}];

        % Loop over categories
        for ic = 1:ncats
            
            if ic == 1
                % if first, include stats
                outcome_cell = [outcome_cell;{cat_labels{ic}},...
                    {sprintf('%d (%1.2f%%)',tab2(2,ic),prct_tab2(2,ic)*100)},...
                    {sprintf('%d (%1.2f%%)',tab2(1,ic),prct_tab2(1,ic)*100)},...
                    {sprintf('%1.1f (%1.1f-%1.1f)',or,ci(1),ci(2))},...
                    {sprintf('%s',simple_p_text(p))}];
                
            else
                outcome_cell = [outcome_cell;{cat_labels{ic}},...
                    {sprintf('%d (%1.2f%%)',tab2(2,ic),prct_tab2(2,ic)*100)},...
                    {sprintf('%d (%1.2f%%)',tab2(1,ic),prct_tab2(1,ic)*100)},...
                    {''},...
                    {''}];
                
            end
            
            
            
        end
        
        
        
    end
    
    
    
end

T = cell2table(outcome_cell,'VariableNames',{'Predictor','Goodoutcome','Pooroutcome','Statistic','P-value'})


end