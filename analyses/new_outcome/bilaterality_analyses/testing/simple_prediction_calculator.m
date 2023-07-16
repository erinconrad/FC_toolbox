function [preds,scores] = simple_prediction_calculator(features,classNames,coefs)

% Get score
scores = 1./(1+exp(-(features*coefs(2)+coefs(1))));

preds = cell(length(features),1);
preds(scores>=0.5) = classNames(2);
preds(scores<0.5) = classNames(1);

end