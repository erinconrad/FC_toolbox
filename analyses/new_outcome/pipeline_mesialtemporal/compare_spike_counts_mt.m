function compare_spike_counts_mt

montage = 1;
montage_mt = 3;
which_sleep_stage = 1;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
subplot_path = [plot_folder,'ai_subplots/'];
if ~exist(subplot_path,'dir')
    mkdir(subplot_path)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;

mt_data = load([inter_folder,'mt_out.mat']);
mt_data = mt_data.out;
names = data.all_names;
assert(isequal(names,mt_data.all_names))

%% Get spikes
mt_spikes = mt_data.all_spikes(:,montage_mt,which_sleep_stage);
spikes = data.all_spikes(:,montage,which_sleep_stage);

%% Get labels
labels = data.all_labels(:,1);
mt_labels = mt_data.all_labels(:,2);

%% Fix mt labels
mt_labels = cellfun(@(X) strrep(X,'-CAR',''),mt_labels,'UniformOutput',false);
mt_labels(cellfun(@isempty,mt_labels)) = {num2cell([])};


%% Find the labels that match mt labels
match = cellfun(@(x,y) ismember(x,y),labels,mt_labels,'UniformOutput',false);
reduced_labels = cellfun(@(x,y) x(y),labels,match,'UniformOutput',false);
reduced_labels(cellfun(@isempty,reduced_labels)) = {num2cell([])};% dumb thing

% Make sure the labels perfectly match
assert(all(cellfun(@(x,y) isequal(x,y),mt_labels,reduced_labels)))

% now reduce the spikes
reduced_spikes = cellfun(@(x,y) x(y),spikes,match,'UniformOutput',false);

if 1
    figure
   % d = find(strcmp(names,'HUP225'));
    for i = 1:length(mt_labels)
        %table(reduced_labels{i},reduced_spikes{i},mt_labels{i},mt_spikes{i})
        if isempty(mt_spikes{i}),continue;end
        r = corr(reduced_spikes{i},mt_spikes{i},'rows','pairwise');
        plot(reduced_spikes{i},mt_spikes{i},'o')
        hold on
        text(reduced_spikes{i},mt_spikes{i},reduced_labels{i})
        xlim = [min([reduced_spikes{i};mt_spikes{i}]) max([reduced_spikes{i};mt_spikes{i}])];
        ylim = [min([reduced_spikes{i};mt_spikes{i}]) max([reduced_spikes{i};mt_spikes{i}])];
        plot(xlim,ylim,'k--')
        title(sprintf('%s r = %1.2f',names{i},r))
        xlabel('Old')
        ylabel('New')
        pause
        hold off
    end
end

end