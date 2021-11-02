function avg_across_times

%{
Pseudocode:
1) Get spikes
2) Remove spikes in seizures
3) average all spike counts together across time
%}

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/avg_spikes_over_time/'];
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

%% Loop over patients
for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get main things
    name = summ.name;
    spikes = summ.spikes;
    labels = summ.labels;
    soz_labels = summ.soz.labels;
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    spikes = spikes(~ekg,:); % spike rate #spikes/elec/one minute block (spikes/elec/min)
    labels = labels(~ekg);
    
    %% Average spikes over all times
    avg_spikes = nanmean(spikes,2); % this will be spikes/elec/min
    
    %% prep out file
    out(p).name = name;
    out(p).avg_spikes = avg_spikes;
    out(p).labels = labels;
    out(p).soz_labels = soz_labels;
    
    if 0
        [sorted_spikes,I] = sort(avg_spikes,'descend');
        table(labels(I),sorted_spikes)
        pause
    end
    
end

%% Save out file
save([out_folder,'out.mat'],'out');


end