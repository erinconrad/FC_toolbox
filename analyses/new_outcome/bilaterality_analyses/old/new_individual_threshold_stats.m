function out = new_individual_threshold_stats(scores,true_labels,positive_label,desired_threshold)

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

% Only good for 2 classes
assert(length(unique(true_labels))==2)
unique_labels = unique(true_labels);

pred_labels = cell(length(true_labels),1);
pred_labels(scores > desired_threshold) = {positive_label};
pred_labels(scores <= desired_threshold) = (unique_labels(~ismember(unique_labels,positive_label)));

out = confusion_matrix(pred_labels,true_labels,0);


end