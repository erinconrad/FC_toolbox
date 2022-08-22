function find_atlas_laterality_disagreement

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;

 %% get variables of interest
aal_atlas_names = data.aal_names;
brainnetome_atlas_names = data.brainnetome_names;
aal = data.all_aal;
labels = data.all_labels;
anatomy = data.all_anatomy;
brainnetome = data.all_brainnetome;
npts = length(aal);


%% Prep cell containing mismatches
mismatches = {};
for ip = 1:npts
    % Get current atlas regions
    curr_aal = aal{ip};
    curr_brainnetome = brainnetome{ip};
    curr_anatomy = anatomy{ip};
    curr_labels = labels{ip};
    
    % Lateralize these
    aal_lats = lateralize_regions_simple(curr_aal);
    brainnetome_lats = lateralize_regions_simple(curr_brainnetome);
    
    any_empty = any([cellfun(@isempty,aal_lats),cellfun(@isempty,brainnetome_lats)],2);
    aal_lats(any_empty) = [];
    brainnetome_lats(any_empty) = [];
    curr_anatomy(any_empty) = [];
    curr_labels(any_empty) = [];
      
    % Look for disageements
    disagree = find(cellfun(@(x,y) ~strcmp(x,y),aal_lats,brainnetome_lats));
    
    % get information about the disagreement
    for id = 1:length(disagree)
        cd = disagree(id);
        mismatches = [mismatches;ip,curr_labels(cd),curr_anatomy(cd),aal_lats(cd),brainnetome_lats(cd)];
    end
end

T = cell2table(mismatches,'VariableNames',{'Pt','label','anatomy','AAL','brainnetome'})

end