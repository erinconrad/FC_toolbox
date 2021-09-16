function [p,post_hoc_p,which_groups] = non_para_circ_stats(data)

data = data';
data(any(isnan(data),2),:) = []; % remove rows with any nans
ngroups = size(data,2);
npts = size(data,1);
post_hoc_p = nan(nchoosek(ngroups,2),1);
which_groups = nan(nchoosek(ngroups,2),2);

if npts == 2
    p = signrank(data);
    
elseif npts > 2
    
    p = friedman(data,1,'off');
    
    % post-hoc tests
    if p < 0.05
        pcount = 0;
        
        % do signrank for each subgroup
        for j = 1:ngroups
            for k = j+1:ngroups
                pcount = pcount+1;
                sp = signrank(data(:,j),data(:,k));
                post_hoc_p(pcount) = sp;
                which_groups(pcount,:) = [j k];
            end
        end
        
    end
    
else
    error('what');
end

end