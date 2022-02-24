function auc = predict_soz_fc(T)

prop_train = 0.67;

%% Divide into training and testing set
training = randsample(length(T.pt_vec),floor(length(T.pt_vec)*prop_train));
training_idx = ismember(T.pt_vec,training);
testing_idx = ~ismember(T.pt_vec,training);

T_train = T(training_idx,:);
T_test = T(testing_idx,:);

% Confirm that I separated by patients
train_pts = unique(T_train.pt_vec);
test_pts = unique(T_test.pt_vec);
assert(isempty(intersect(train_pts,test_pts)))

%% Train model with glm model
%
glm = fitglm(T_train,'soz_vec ~ spike_vec + ns_vec',...
'Distribution','Poisson','Link','log');
%}
%{
glm = fitglm(T_train,'soz_vec ~ spike_vec',...
'Distribution','Poisson','Link','log');
%}

%% Test model on testing data
params = glm.CoefficientNames;
asum = zeros(size(T_test,1),1);
asum = asum + glm.Coefficients.Estimate(1);
for p = 2:length(params)
    est = glm.Coefficients.Estimate(p);
    asum = asum +  T_test.(params{p})*est;
end

classification = logistic(asum);

all_soz = T_test.soz_vec==1;
all_no_soz = T_test.soz_vec==0;
class_soz = classification(all_soz);
class_no_soz = classification(all_no_soz);
[roc,auc,disc,disc_I] = calculate_roc(class_no_soz,class_soz,1e3);


end

function out = logistic(x)

out = 1./(1+exp(-x));

end
