function [p,post_hoc_p,which_groups] = non_parametric_stats(data)

ngroups = size(data,1);
npts = size(data,2);
post_hoc_p = nan(nchoosek(ngroups,2),1);
which_groups = nan(nchoosek(ngroups,2),2);

if npts == 2
    p = signrank(data');
    
elseif npts > 2
    p = skillmack(data);
    
    % post-hoc tests
    if p < 0.05
        pcount = 0;
        
        % do paired t-tests for each subgroup
        for j = 1:ngroups
            for k = 1:j-1
                pcount = pcount+1;
                sp = signrank(data(j,:)',data(k,:)');
                post_hoc_p(pcount) = sp;
                which_groups(pcount,:) = [j k];
            end
        end
        
    end
    
    
    
else
    error('what');
end


end