function [npv,ppv,npred,totaln,mat,desired_threshold,acc,sens,spec] = individual_threshold_stats(scores,true_labels,desired_threshold)

if isempty(scores)
    npv = nan;
    ppv = nan;
    npred = nan;
    totaln = nan;
    mat = [nan,nan;nan,nan];
    desired_threshold = nan;
    acc = nan;
    sens = nan;
    spec = nan;
    return
end

if isempty(desired_threshold)
    
    %{
    nsoz = sum(true_labels == 1);
    
    % Calculate the number of predicted SOZ for each possible threshold in T
    n_pred_soz = nan(length(T),1);
    for it = 1:length(T)
        n_pred_soz(it) = sum(scores > T(it));
    end

    % Find the index that minimizes the difference between predicted number of
    % SOZ and true number
    [~,I] = min(abs(nsoz-n_pred_soz));

    % Get corresponding threshold
    desired_threshold = T(I);
    %}
end

if isnumeric(true_labels)
    true_labels = arrayfun(@num2str,true_labels,'uniformoutput',false);
end

pred_labels = cell(length(true_labels),1);
pred_labels(scores > desired_threshold) = {'1'};
pred_labels(scores <= desired_threshold) = {'0'};

out = confusion_matrix(pred_labels,true_labels,0);
ppv = out.ppv;
npv = out.npv;
mat = out.mat;
acc = out.accuracy;
npred = sum(scores > desired_threshold);
totaln = length(scores);
sens = out.sensitivity;
spec = out.specificity;

end