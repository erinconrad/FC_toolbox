function density = estimate_coverage_density(locs,r,labels)

%% Overview
%{
This function takes a set of points and calculates the density at each one.
To do this, I follow the algorithm from ArcGIS:
https://doc.arcgis.com/en/insights/latest/analyze/calculate-density.htm
which cites Silverman (1986, p. 76, Equation 4.5).

I modify this to ignore the electrode at the point of interest (because
otherwise density would be zero). 

parameter r is search radius

%}


%% Parameters
W = 1; % weight (allows you to weight different point differently)

%% Initialize output variable
nelecs = size(locs,1);
density = nan(nelecs,1);

%% Fix for weird locs
locs(any(locs > 1e5,2),:) = nan;

%% Loop over all elecs to calculate density
for in = 1:nelecs
    
    %% Initialize density
    curr_dens = 0;
    
    %% Loop over the remaining elecs
    for jn = 1:nelecs
        
        % skip itself
        if in == jn, continue; end
        if any(isnan(locs(jn,:))), continue; end
          
        % get distance between i and j elecs
        %d = sqrt(sum((locs(in,:)-locs(jn,:)).^2));
        d = vecnorm(locs(in,:)-locs(jn,:));
        
        % skip points outside the search radius
        if d > r, continue; end
        
        % add kernel function for current elec
        curr_dens = curr_dens + 3/pi * W * (1 - (d/r)^2)^2;
        
    end
    
    %% Finalize density calculation
    curr_dens = curr_dens * 1/r^2;
    
    %% Output
    density(in) = curr_dens;
    
end

% visualize
if 0
    figure
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,density,'filled',...
        'markeredgecolor','k','linewidth',2);
    hold on
    text(locs(:,1),locs(:,2),locs(:,3),labels,...
        'horizontalalignment','center')
    c = colorbar;
    xlabel('X position (mm)')
    ylabel('Y position (mm)')
    zlabel('Z position (mm)')
    ylabel(c,'Density (1/mm^2)')
    set(gca,'fontsize',15)
    title('Electrode coverage density')
end


end