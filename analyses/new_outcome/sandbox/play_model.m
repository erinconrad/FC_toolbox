function play_model

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [results_folder,'analysis/new_outcome/data/'];
out_folder = [results_folder,'analysis/new_outcome/out/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Load data file
data = load([data_folder,'main_out.mat']);
data = data.out;

%% Get variables
ilae = data.all_two_year_ilae;
npts = length(ilae);
spikes = data.all_spikes;
fc = data.all_fc;
stereo = data.all_stereo;

%% bin outcomes into good and bad
no_outcome = cellfun(@isempty,ilae);
good_outcome = cellfun(@(x) contains(x,'1'),ilae);
bad_outcome = cellfun(@(x) any(contains(x,{'2','3','4','5','6'})),ilae);
outcome = nan(npts,1);
outcome(good_outcome) = 1;
outcome(bad_outcome) = 0;

% Make sure I have captured everyone
assert(sum(no_outcome)+sum(good_outcome)+sum(bad_outcome) == npts);

%% NUmber of electrodes
nelecs = cellfun(@length,spikes);

%% Mean spike rates (note this includes bad spikes!!!)
avg_spikes = cellfun(@nanmean,spikes);

%% Mean fc
avg_fc = cellfun(@(x) nanmean(x,'all'),fc);

%% Does mean FC vary by implant strategy?
% yes
if 0
figure
unpaired_plot(avg_fc(stereo==1),avg_fc(stereo==0),{'Stereo','G/S/D'},'Average correlation')
end

%% Prepare structure of predictors
count = 1;
% Number of electrodes
predictors(count).name = 'Number of electrodes';
predictors(count).type = 'continuous';
predictors(count).X = nelecs;
count = count+1;

% Stereo
stereo_cell = cell(npts,1);
stereo_cell(stereo==1) = {'Stereo'};
stereo_cell(~stereo) = {'G/S/D'};
predictors(count).name = 'Implant type';
predictors(count).type = 'binary';
predictors(count).X = stereo_cell;
count = count+1;

% spikes
predictors(count).name = 'Spikes';
predictors(count).type = 'continuous';
predictors(count).X = avg_spikes;
count = count+1;

% FC
predictors(count).name = 'Connectivity';
predictors(count).type = 'continuous';
predictors(count).X = avg_fc;
count = count + 1;


%% Table of univariate comparisons
T = outcome_table(outcome,predictors);
writetable(T,[out_folder,'outcome_table.csv'])

end