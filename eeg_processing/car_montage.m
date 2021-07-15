function values = car_montage(values,which_chs)

% Do car (just average the non-skip chs)
values = values - repmat(nanmean(values(:,which_chs),2),1,size(values,2));


end