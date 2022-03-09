function out = sleep_conf_matrix(class_soz,class_no_soz,disc)

tp = sum(class_soz>disc);
fn = sum(class_soz<disc);
fp = sum(class_no_soz>disc);
tn = sum(class_no_soz<disc);

sensitivity = tp/(tp+fn); % what % of those truly positive test positive?
specificity = tn/(tn+fp); % what % of those truly negative test negative?
ppv = tp/(tp+fp); % what % of those predicted to be positive are positive?
npv = tn/(tn+fn); % what % of those predicted to be negative are negative?

mat = [tp fn;fp tn];

out.mat = mat;
out.sensitivity = sensitivity;
out.specificity = specificity;
out.ppv = ppv;
out.npv = npv;

end