function combine_subplots

%% Parameters
which_things = {'spikes_1','pearson_1'};

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
subplot_path = [plot_folder,'ai_subplots/'];
out_path = [plot_folder,'ai_plots/'];
if ~exist(out_path,'dir')
    mkdir(out_path)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;
names = data.all_names;
npts = length(names);

%% Loop over patients
for ip = 1:npts
    name = names{ip};

    % Prep out figure
    mf = figure;
    set(gcf,'position',[10 10 1000 800]);
    t = tiledlayout(length(which_things),1);
    for it = 1:length(which_things)
        % try to load the relevant subfigure
        subfig = [subplot_path,name,'_',which_things{it},'.fig'];
        if exist(subfig,'file') ~=0
            f = openfig(subfig);
            ax2 = f.Children;
            ax2.Parent=t;
            ax2.Layout.Tile = it;
        end
        

    end
    print(mf,[out_path,name],'-dpng')
    close all
end

end