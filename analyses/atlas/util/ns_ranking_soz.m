function ns_ranking_soz(z_score_ns,sozs,atlas_nums,atlas_names,pt_names)

npts = length(sozs);
nums = 1:length(atlas_nums);

all_rankings = nan(npts,2);

if 0
    turn_nans_gray(z_score_ns)
    hold on
    for ip = 1:npts
        soz = sozs{ip};
        soz_indices = (ismember(atlas_nums,soz));
        plot(repmat(ip,sum(soz_indices),1),find(soz_indices),'rx','linewidth',2)
    end
    yticks(1:size(z_score_ns,1))
    yticklabels(atlas_names)
    
end

% Loop over patients
for ip = 1:npts
        
    % Add the soz locs
    soz = sozs{ip};
    soz_indices = (ismember(atlas_nums,soz));
    
    % remove nan zscores
    curr_ns = z_score_ns(:,ip);
    nan_idx = isnan(curr_ns);
    
    soz_indices(nan_idx) = [];
    curr_ns(nan_idx) = [];
        
    % order the z_score_ns
    [~,I] = sort(curr_ns,'descend');
    
    % convert to rankings
    ranks = 1:length(curr_ns);
    ranks(I) = ranks;
    
    % Get ranks of SOZ elecs
    soz_ranks = nanmedian(ranks(soz_indices));
    all_rankings(ip,:) = [soz_ranks,nanmedian(ranks)];

    % Plot
    if 0
        plot(nums,z_score_ns(:,ip),'o','linewidth',2)
        hold on
        plot(nums(soz_indices),z_score_ns(soz_indices,ip),'ro','linewidth',2);
        title(pt_names{ip})
        pause
        hold off
    end
end

figure
max_pos = max(all_rankings(:));
min_pos = min(all_rankings(:));
plot(all_rankings(:,1),all_rankings(:,2),'o','linewidth',2);
hold on
plot([min_pos max_pos],[min_pos max_pos],'k--','linewidth',2);



end