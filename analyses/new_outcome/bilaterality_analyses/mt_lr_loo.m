function T = mt_lr_loo(T,features)

% seed rng
rng(0)

do_plot = 1;
method = 'bag';
ncycles = 1e3;
response = 'soz_lats';%'soz_lats';%'outcome';
pca_perc = 70;
which_outcome = 'engel';
which_year = 1;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];

%% Remove patients without response
empty_class = cellfun(@isempty,T.(response));
T(empty_class,:) = [];
npts = size(T,1);

% initialize confusion matrix
classes = unique(T.(response));
nclasses = length(classes);
C = zeros(nclasses,nclasses); % left, right, bilateral
all_pred = cell(npts,1);
Cnull = zeros(nclasses,nclasses);

%% Do leave-one-patient-out classifier to predict laterality
for i = 1:npts
    
    % split into training and testing
    Ttrain = T([1:i-1,i+1:end],:); % training all but current
    Ttest = T(i,:); % testing current

    % make sure they're distinct
    assert(isempty(intersect(Ttrain.names,Ttest.names)))
    features_minus_null = features;
    features_minus_null(strcmp(features_minus_null,'nelecs_1_car')) = [];

    % train classifier    
    tc = lr_classifier(Ttrain,method,features_minus_null,response,pca_perc,ncycles);
    tcnull = lr_classifier(Ttrain,method,{'spikes_1_car'},response,pca_perc,ncycles);

    % make prediction on left out
    pred = tc.predictFcn(Ttest);
    all_pred{i} = pred{1};
    prednull = tcnull.predictFcn(Ttest);
    
    % compare to true
    true = Ttest.(response);

    % which row to add to confusion matrix (the true value)
    which_row = find(strcmp(true,classes));

    % which column to add to confusion matrix (the predicted value)
    which_column = find(strcmp(pred,classes));

    C(which_row,which_column) = C(which_row,which_column) + 1;

    which_column_null = find(strcmp(prednull,classes));
    Cnull(which_row,which_column_null) = Cnull(which_row,which_column_null) + 1;
    
    
end

%% Calculate accuracy and balanced accuracy
accuracy = sum(diag(C))/sum(C(:));
accuracy_null = sum(diag(Cnull))/sum(Cnull(:));

% Balanced accuracy is the average across all classes of the number of 
% data accurately predicted belonging to class m divided by the number of
% data belonging to class m
recall = nan(nclasses,1);
recall_null = nan(nclasses,1);
for i = 1:nclasses
    tp = C(i,i);
    fn = sum(C(i,~ismember(1:nclasses,i))); 
    recall(i) = tp/(tp+fn); % tp is number correctly predicted to be in class, tp + fn is everyone in the class

    tpnull = Cnull(i,i);
    fnnull = sum(Cnull(i,~ismember(1:nclasses,i)));
    recall_null(i) = tpnull/(tpnull+fnnull);
end
balanced_accuracy = mean(recall);
balanced_accuracy_null = mean(recall_null);

%% Double check accuracy another way
assert(sum(cellfun(@(x,y) strcmp(x,y),all_pred,T.(response)))/length(all_pred)==accuracy)
if 0
    table(T.(response),all_pred)
end

%% Plot
figure
set(gcf,'position',[10 10 1400 500])
tiledlayout(1,2,"TileSpacing",'tight','padding','tight')

%% confusion matrix
% Map numbers onto 0 to 1
nexttile
new_numbers = map_numbers_onto_range(C,[1 0]);
Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
D = diag(new_numbers);
Dcolor = [repmat(D,1,2),ones(length(D),1)];
Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
imagesc(Ccolor)

xticks(1:nclasses)
xticklabels((classes))
yticks(1:nclasses)
yticklabels((classes))
xlabel('Predicted')
ylabel('True')
hold on
for i = 1:nclasses
    for j = 1:nclasses
        text(i,j,sprintf('%d',C(j,i)),'horizontalalignment','center','fontsize',20)
    end
end

title(sprintf('Accuracy: %1.1f%%\nBalanced accuracy: %1.1f%%',...
    accuracy*100,balanced_accuracy*100))
set(gca,'fontsize',20)

nexttile
new_numbers = map_numbers_onto_range(Cnull,[1 0]);
Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
D = diag(new_numbers);
Dcolor = [repmat(D,1,2),ones(length(D),1)];
Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
imagesc(Ccolor)

xticks(1:nclasses)
xticklabels((classes))
yticks(1:nclasses)
yticklabels((classes))
xlabel('Predicted')
ylabel('True')
hold on
for i = 1:nclasses
    for j = 1:nclasses
        text(i,j,sprintf('%d',Cnull(j,i)),'horizontalalignment','center','fontsize',20)
    end
end

title(sprintf('Null Accuracy: %1.1f%%\nBalanced accuracy: %1.1f%%',...
    accuracy_null*100,balanced_accuracy_null*100))
set(gca,'fontsize',20)

print(gcf,[plot_folder,'model'],'-dpng')

%% Save data
T = addvars(T,all_pred,'NewVariableNames','pred_lat','After','surg_lat');
writetable(T,[plot_folder,'model_pred.csv'])

%{
%% Now predict outcome

% Remove patients with missing data
no_outcome = cellfun(@isempty,T.(response));
no_surg = ~strcmp(T.surgery,'Laser ablation') & ~strcmp(T.surgery,'Resection');
not_temporal = ~strcmp(T.surg_loc,'temporal');

remove = no_outcome | no_surg | not_temporal;
oT = T;
oT(remove,:) = [];

% Parse actual outcome
outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
outcome = cellfun(@(x) parse_outcome_new(x,which_outcome),oT.(outcome_name),'UniformOutput',false);


%{\
% Predict good outcome if predicted to be unilateral and agrees with surg
% laterality
%unilateral_pred = strcmp(oT.pred_lat,'left') | strcmp(oT.pred_lat,'right');
agree = cellfun(@(x,y) strcmp(x,y),oT.surg_lat,oT.pred_lat);
%agree = cellfun(@(x,y) strcmp(x,y),oT.soz_lats,oT.pred_lat);
%pred_good = unilateral_pred & agree; 
pred_good = agree;

% Predict bad outcome if predicted to be bilateral OR if disagrees with
% surg laterality
bilateral_pred = strcmp(oT.pred_lat,'bilateral');
disagree = ~agree;
%pred_bad = bilateral_pred | disagree;
pred_bad = disagree;

assert(isequal(pred_good,~pred_bad))

%% Confusion matrix for predicted and true outcome
% Make a confusion matrix for outcome
Cout(1,1) = sum(strcmp(outcome,'good') & pred_good  == 1);
Cout(1,2) = sum(strcmp(outcome,'good') & pred_good == 0);
Cout(2,1) = sum(strcmp(outcome,'bad') & pred_good  == 1);
Cout(2,2) = sum(strcmp(outcome,'bad') & pred_good == 0);
accuracy_out = sum(diag(Cout))/sum(Cout(:));




if do_plot
    figure
    set(gcf,'position',[10 10 1000 400])
    tiledlayout(1,2)

    %% confusion matrix
    nexttile
    % Map numbers onto 0 to 1
    new_numbers = map_numbers_onto_range(C,[1 0]);
    Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
    D = diag(new_numbers);
    Dcolor = [repmat(D,1,2),ones(length(D),1)];
    Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
    imagesc(Ccolor)
    
    xticks(1:nclasses)
    xticklabels((classes))
    yticks(1:nclasses)
    yticklabels((classes))
    xlabel('Predicted')
    ylabel('True')
    hold on
    for i = 1:nclasses
        for j = 1:nclasses
            text(i,j,sprintf('%d',C(j,i)),'horizontalalignment','center','fontsize',20)
        end
    end
    
    title(sprintf('Accuracy: %1.1f%%\nBalanced accuracy: %1.1f%%',...
        accuracy*100,balanced_accuracy*100))
    set(gca,'fontsize',20)

    %% outcome for those with agreement vs those without
    nexttile
    % Map numbers onto 0 to 1
    new_numbers = map_numbers_onto_range(Cout,[1 0]);
    nclasses = 2;
    classes = {'Good','Poor'};
    Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
    D = diag(new_numbers);
    Dcolor = [repmat(D,1,2),ones(length(D),1)];
    Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
    imagesc(Ccolor)
    
    xticks(1:nclasses)
    xticklabels((classes))
    yticks(1:nclasses)
    yticklabels((classes))
    xlabel('Predicted')
    ylabel('True')
    hold on
    for i = 1:nclasses
        for j = 1:nclasses
            text(i,j,sprintf('%d',Cout(j,i)),'horizontalalignment','center','fontsize',20)
        end
    end
    
    title(sprintf('Accuracy: %1.1f%%',accuracy_out*100))
    set(gca,'fontsize',20)

    %% Null accuracy

end
%}

end