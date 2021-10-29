function [p,post_hoc_p,which_groups] = non_parametric_stats(data,p_or_unp)

ngroups = size(data,1);
npts = size(data,2);
post_hoc_p = nan(nchoosek(ngroups,2),1);
which_groups = nan(nchoosek(ngroups,2),2);

switch p_or_unp
    case 'paired'
        if ngroups == 2
            p = signrank(data');

        elseif ngroups > 2
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

        
    case 'unpaired'
        if ngroups == 2
            p = ranksum(data(1,:)',data(2,:)');
        elseif ngroups > 2
            p = kruskalwallis(data');
            
            % post-hoc tests
            if p < 0.05
                pcount = 0;

                % do paired t-tests for each subgroup
                for j = 1:ngroups
                    for k = 1:j-1
                        pcount = pcount+1;
                        sp = ranksum(data(j,:)',data(k,:)');
                        post_hoc_p(pcount) = sp;
                        which_groups(pcount,:) = [j k];
                    end
                end

            end
        end
        %}
end