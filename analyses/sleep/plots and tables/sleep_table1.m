function T = sleep_table1

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

names = cell(npts,1);
age_onset = nan(npts,1);
age_implant = nan(npts,1);
sex = cell(npts,1);
nelecs = nan(npts,1);
any_grids = nan(npts,1);
rate = nan(npts,1);
duration = nan(npts,1);
lat = cell(npts,1);
loc = cell(npts,1);

%% Loop over files
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Basic demographics
    names{p} = summ.name;
    clinical = summ.clinical;
    age_onset(p) = clinical.age_onset;
    sex{p} = clinical.sex;
    age_implant(p) = clinical.age_implant;
    
    %% Electrodes
    labels = summ.labels;
    ekg = find_non_intracranial(labels);
    labels = labels(~ekg);
    
    nelecs(p) = length(labels);
    any_grids(p) = decide_if_any_grids_or_strips(labels);
    
    %% Spike rate
    spikes = summ.spikes;
    rate(p) = nanmean(spikes,'all');
    
    
    %% Duration
    duration(p) = summ.times(end)/3600/24;
    
    %% SOZ localization
    lat{p} = summ.soz.lat;
    loc{p} = summ.soz.loc;
    
    
end

T = table(names,sex,age_onset,age_implant,nelecs,any_grids,...
    duration,rate,loc,lat);


end


function any_grid = decide_if_any_grids_or_strips(labels)

any_grid = 0;

for i = 1:length(labels)
    curr = labels{i};
    if contains(curr,'G')
        B = regexp(curr,'\d*','Match');
        if isempty(B), continue; end
        B = str2num(B{1});
        if B > 12 % if it's *G13 or higher then it's a grid
            any_grid = 1;
            break
        end
    end
end

end