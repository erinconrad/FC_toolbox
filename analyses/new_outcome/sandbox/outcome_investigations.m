function outcome_investigations(T,features)

which_outcome = 'engel';
which_year = 1;
rm_non_temporal = 1;

%% Remove non temporal patients
if rm_non_temporal
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];
end


%% Parse good vs bad outcome
% Parse actual outcome
outcome_name = [which_outcome,'_yr',sprintf('%d',which_year)];
outcome = cellfun(@(x) parse_outcome_new(x,which_outcome),T.(outcome_name),'UniformOutput',false);
T.outcome = outcome;

%% Compare outcomes by laterality
if 1
figure
heatmap(T,'soz_lats',outcome_name)
end

%% Compare outcomes by spike AI
if 0
feature = 'spikes car sleep';
figure
boxplot_with_points(T.(feature),T.(outcome_name),0)
end


end