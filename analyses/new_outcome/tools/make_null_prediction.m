function [N,accuracy,balanced_accuracy] = make_null_prediction(C)

totalN = sum(C,'all');
true_probs = sum(C,2)/totalN;
N = zeros(size(C,1),size(C,2));
for r = 1:size(C,1)
    nrow = sum(C(r,:)); % keep the number in the row the same, but I am going to shuffle the column
    for j = 1:nrow % loop over number of things in the row
        b = weighted_rand(size(C,2),true_probs); % pick a random column, probability weighted by true probabilities
        N(r,b) = N(r,b) + 1; % put it in that column and the original row
    end
end

accuracy = sum(diag(N))/sum(N(:));
% Balanced accuracy is the average across all classes of the number of 
% data accurately predicted belonging to class m divided by the number of
% data belonging to class m
recall = nan(size(C,1),1);
for i = 1:size(C,1)
    tp = C(i,i);
    fn = sum(C(i,~ismember(1:size(C,1),i))); 
    recall(i) = tp/(tp+fn); % tp is number correctly predicted to be in class, tp + fn is everyone in the class
end
balanced_accuracy = mean(recall);


end

function b = weighted_rand(a,weights)

b = randsample(a,1,true,weights);

end