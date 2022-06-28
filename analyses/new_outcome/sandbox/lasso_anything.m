function c = lasso_anything(yBinom,X)

%{
Should consider if the threshold of 0.5 makes sense
%}

%% Separate into training and testing data
rng default    % Set the seed for reproducibility
c = cvpartition(yBinom,'HoldOut',0.3);
idxTrain = training(c,1);
idxTest = ~idxTrain;
XTrain = X(idxTrain,:);
yTrain = yBinom(idxTrain);
XTest = X(idxTest,:);
yTest = yBinom(idxTest);


%% Perform lasso regularization with logistic regression and 3-fold cross-validation on training data
[B,FitInfo] = lassoglm(XTrain,yTrain,'binomial','CV',3);
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];


%% Test on testing data
yhat = glmval(coef,XTest,'logit');

%% Pick threshold of 0.5 and generate confusion chart
yhatBinom = (yhat>=0.5);
c = confusionmat(yTest,yhatBinom);


end