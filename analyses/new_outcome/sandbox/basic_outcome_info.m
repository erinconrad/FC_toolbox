function basic_outcome_info

use_raw = 0;
which_year = 2;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

if use_raw
    
    %% Load data file
    pt = load([data_folder,'pt.mat']);
    pt = pt.pt;
    npts = length(pt);
    
    ilae = cell(npts,2);
    engel = cell(npts,2);
    
    for ip = 1:npts
        ilae(ip,:) = pt(ip).clinical.ilae';
        engel(ip,:) = pt(ip).clinical.engel';
        
        assert(isequal(pt(ip).clinical.engel_years,[1 2]))
        assert(isequal(pt(ip).clinical.ilae_years,[1 2]))
    end
    
    %% Just make two year
    ilae = ilae(:,which_year);
    engel = engel(:,which_year);
    
    
    %% Get all unique outcomes in both sets
    [engel_cats,~,ic_engel] = unique(engel);
    [ilae_cats,~,ic_ilae] = unique(ilae);

    %% Make a summary table
    engel_counts = cellfun(@(x) sum(strcmp(x,engel)),engel_cats);
    ilae_counts = cellfun(@(x) sum(strcmp(x,ilae)),ilae_cats);
    
    fprintf('\nThere are %d patients with %d year ilae data.\n',sum(ilae_counts(2:end)),which_year);
    figure
    nexttile
    bar(ilae_counts(2:end))
    xticklabels(ilae_cats(2:end))
    
    nexttile
    bar(engel_counts(2:end))
    xticklabels(engel_cats(2:end))
    
    
else

    %% Load data file
    data = load([inter_folder,'main_out.mat']);
    data = data.out;

    %% outcomes
    ilae = data.all_two_year_ilae;
    engel = data.all_two_year_engel;

    %% Get all unique outcomes in both sets
    [engel_cats,~,ic_engel] = unique(engel);
    [ilae_cats,~,ic_ilae] = unique(ilae);

    %% Make a summary table
    engel_counts = cellfun(@(x) sum(strcmp(x,engel)),engel_cats);
    ilae_counts = cellfun(@(x) sum(strcmp(x,ilae)),ilae_cats);

    table(engel_cats,engel_counts)
    table(ilae_cats,ilae_counts)
    
    
    
end

end