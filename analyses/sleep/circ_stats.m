function [p,stats,post_hoc_p,which_groups] = circ_stats(data)

ngroups = size(data,1);
npts = size(data,2);
post_hoc_p = nan(nchoosek(ngroups,2),1);
which_groups = nan(nchoosek(ngroups,2),2);

if npts == 2
    [~,p,~,stats] = ttest(data');
    
elseif npts > 2
    t = [(1:size(data,2))',data'];
    t(any(isnan(t),2),:) = [];
    T = array2table(t);
    if size(T,2) == 5
        rm = fitrm(T,'t2-t5~t1');
    elseif size(T,2) == 3
        rm = fitrm(T,'t2-t3~t1');
    end
    ranovatbl = ranova(rm);
    p = ranovatbl.pValue(1);
    stats = ranovatbl.F(1);
    
    % post-hoc tests
    if p < 0.05
        pcount = 0;
        
        % do paired t-tests for each subgroup
        for j = 1:ngroups
            for k = 1:j-1
                pcount = pcount+1;
                [~,sp] = ttest(data(j,:)',data(k,:)');
                post_hoc_p(pcount) = sp;
                which_groups(pcount,:) = [j k];
            end
        end
        
    end
    
    
    % test that when I just take two measurements this gives basically the
    % same result as a paired t test
    if 0
        rm2 = fitrm(T,'t2-t3~t1');
        ranovatbl2 = ranova(rm2);
        alt_data = t(:,[2 3]);
        [~,ptest] = ttest(alt_data(:,1),alt_data(:,2));
        
        
    end
    
    
else
    error('what');
end

end