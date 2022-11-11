function all_pred = mt_lr_loo(T,features)

do_plot = 1;

%% Remove non temporal?
non_temporal = ~strcmp(T.soz_locs,'temporal');
%T(non_temporal,:) = [];

%% Remove bilateral?
bilateral = strcmp(T.soz_lats,'bilateral');
%T(bilateral,:) = [];

npts = size(T,1);

% initialize confusion matrix
classes = unique(T.soz_lats);
nclasses = length(classes);
C = zeros(nclasses,nclasses); % left, right, bilateral

all_pred = cell(npts,1);
all_true = T.soz_lats;
all_outcome = T.outcome;

%% Do classifier to predict laterality
for i = 1:npts
    
    % split into training and testing
    Ttrain = T([1:i-1,i+1:end],:);
    Ttest = T(i,:);

    % make sure they're distinct
    assert(isempty(intersect(Ttrain.names,Ttest.names)))

    % train classifier    
    tc = lt_mr_tree(Ttrain,'bag',features);

    % make prediction on left out
    pred = tc.predictFcn(Ttest);
    all_pred{i} = pred{1};

    % compare to true
    true = Ttest.soz_lats;

    % which row to add to confusion matrix (the true value)
    which_row = find(strcmp(true,classes));

    % which column to add to confusion matrix (the predicted value)
    which_column = find(strcmp(pred,classes));

    C(which_row,which_column) = C(which_row,which_column) + 1;
    
    
end

accuracy = sum(diag(C))/sum(C(:));


%% Now predict outcome
T = addvars(T,all_pred,'NewVariableNames','pred_lat','After','soz_lats');

% Remove patients with missing data
no_outcome = isnan(T.outcome);
no_surg = ~strcmp(T.surgery,'Laser ablation') & ~strcmp(T.surgery,'Resection');
not_temporal = ~strcmp(T.surg_loc,'temporal');

remove = no_outcome | no_surg | not_temporal;
oT = T;
oT(remove,:) = [];

%{\
% Predict good outcome if predicted to be unilateral and agrees with SOZ
% laterality
unilateral_pred = strcmp(oT.pred_lat,'left') | strcmp(oT.pred_lat,'right');
agree = cellfun(@(x,y) strcmp(x,y),oT.soz_lats,oT.pred_lat);
pred_good = unilateral_pred & agree; 

% Predict bad outcome if predicted to be bilateral OR if disagrees with side of SOZ
bilateral_pred = strcmp(oT.pred_lat,'bilateral');
disagree = ~agree;
pred_bad = bilateral_pred | disagree;

assert(isequal(pred_good,~pred_bad))
%}
%{
% Predict good outcome if predict unilateral
pred_good = strcmp(all_pred,'left') | strcmp(all_pred,'right');

% Predict bad outcome if predict bilateral
pred_bad = strcmp(all_pred,'bilateral');
%}


%% Confusion matrix for predicted and true outcome
% Make a confusion matrix for outcome
Cout(1,1) = sum(oT.outcome == 1 & pred_good  == 1);
Cout(1,2) = sum(oT.outcome == 1 & pred_good == 0);
Cout(2,1) = sum(oT.outcome==0 & pred_good  == 1);
Cout(2,2) = sum(oT.outcome==0 & pred_good == 0);
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
    
    title(sprintf('Accuracy: %1.1f%%',accuracy*100))
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

end