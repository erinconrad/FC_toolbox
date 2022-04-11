function sr = calculate_default_search_radius(all_locs)

%{
This function follows arcGIS algorithm to determine default search radius,
based on Silverman's Rule-of-thumb bandwidth
Silverman, B. W. Density Estimation for Statistics and Data Analysis. New York: Chapman and Hall, 1986.
https://desktop.arcgis.com/en/arcmap/latest/tools/spatial-analyst-toolbox/how-kernel-density-works.htm#ESRI_SECTION2_010BE86F10B94294A99F5EF9BF142EB1
%}

npts = length(all_locs);
all_sr = nan(npts,1);

for ip = 1:npts
    locs = all_locs{ip};
    locs(any(locs > 1e5,2),:) = nan; % fix for one weird patient
    
    if isnan(nanmean(locs,'all')), continue; end
    
    % mean center of input points
    mean_center = nanmean(locs,1);
    nelecs = size(locs,1);
    
    if 0
        figure
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k',...
            'markeredgecolor','k','linewidth',2);
        hold on
        scatter3(mean_center(:,1),mean_center(:,2),mean_center(:,3),100,'r','filled');
    end
    
    % get the distance from the mean center for all points
    D = vecnorm(locs - repmat(mean_center,nelecs,1),2,2);
    
    % get median distance
    Dm = nanmedian(D);
    
    % get standard distance
    all_SD = nan(3,1);
    
    % Loop over spatial dimensions
    for d = 1:3
        SD = 0; % set current dimension SD to zero
        for in = 1:nelecs
            if isnan(locs(in,d)), continue; end
            SD = SD + (locs(in,d)-mean_center(d)).^2; % add square distance in that dimesion for that point
        end
        SD = SD/nelecs; % divide by number of elecs
        all_SD(d) = SD;
    end
    SD = sqrt(sum(all_SD)); % add across dimensions, take sqrt
    
    sr = 0.9*min([SD,sqrt(1/log(2))*Dm])*nelecs^(-0.2);
    
    if sr == 0, error('why'); end
    
    all_sr(ip) = sr;
    
end

sr = nanmedian(all_sr);

end