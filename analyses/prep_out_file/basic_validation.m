function basic_validation

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

for p = 1:npts
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get main things
    name = summ.name;
    loc = summ.ana_loc;
    spikes = summ.spikes;
    labels = summ.labels;
    
    %% Get features for soz vs not
    soz = summ.soz.chs;
    chnums = 1:length(labels);
    is_soz = ismember(chnums,soz);
    
    %% Skip subsequent loc analyses if missing
    if sum(cellfun(@(x) isempty(x),loc)) == length(loc) 
        fprintf('\nMissing locs for %s\n',name);
        continue
    end
    
    %% Find and remove non-intracranial
    %{
    MUST REMEMBER TO ADD THIS FOR COA
    %}
    ekg = find_non_intracranial(labels);

    loc = loc(~ekg,:);
    spikes = spikes(~ekg,:);
    labels = labels(~ekg);
    is_soz = is_soz(~ekg);
    
    %% Table showing locs
    fprintf('\nLocations for %s:\n',name);
    table(labels,loc)
    
    %% get top 5 spiking channels
    spikes(isnan(nanmean(spikes,2)),:) = -inf;
    [~,I] = sort(nanmean(spikes,2),'descend');
    fprintf('\nTop 5 spiking channels for %s:\n',name);
    table(labels(I(1:5)),loc(I(1:5)))
    
    %% Shows soz
    fprintf('\nSOZ for %s\n',name);
    table(labels(is_soz),loc(is_soz))
    
    pause
    
    
    
end

end