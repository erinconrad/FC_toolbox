function show_all_outcomes(T)

which_outcomes = {'engel','ilae'};
which_years = [1 2];

figure
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')

%% Remove missing data
no_surg = ~strcmp(T.surgery,'Laser ablation') & ~strcmp(T.surgery,'Resection');
not_temporal = ~strcmp(T.surg_loc,'temporal');
remove = no_surg | not_temporal;
%T(remove,:) = [];

for o = 1:length(which_outcomes)
    for y = which_years

        outcome_name = [which_outcomes{o},'_yr',sprintf('%d',y)];
        outcome = T.(outcome_name);

        % remove empty
        outcome = outcome(cellfun(@(x) ~isempty(x),outcome));
        C = categorical(outcome);
        
        [uni,~,idx] = unique(outcome);
        counts = accumarray(idx(:),1,[],@sum);

        % parse good vs bad
        good_bad = cellfun(@(x) parse_outcome_new(x,which_outcomes{o}),outcome,'UniformOutput',false);
        good = strcmp(good_bad,'good');
        bad = strcmp(good_bad,'bad');

        nexttile
        
        histogram(C(good),uni)
        hold on
        histogram(C(bad),uni)

        title(sprintf('%s year %d',which_outcomes{o},y))
        ylabel('number of patients')
        set(gca,'fontsize',15)

    end
end

end