function [values,car_labels] = car_montage(values,which_chs,labels)

car_labels = labels;

% Do car (just average the non-skip chs)
values = values - repmat(nanmean(values(:,which_chs),2),1,size(values,2));

for i = 1:length(labels)
    car_labels{i} = [labels{i},'-CAR'];
end

end