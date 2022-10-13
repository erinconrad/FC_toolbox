function count_number_with_outcome(which_outcome,which_year)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
new_all_folder = [results_folder,'analysis/intermediate/'];
new_epilepsia_folder = [results_folder,'analysis/intermediate_epilepsia_revision/'];


%% Loop over intermediate files
listing = dir([new_epilepsia_folder,'*.mat']);


%% Prep variables
npts = length(listing);
all_two_year_ilae = cell(npts,1);
all_one_year_ilae = cell(npts,1);
all_two_year_engel = cell(npts,1);
all_one_year_engel = cell(npts,1);
all_surgery = cell(npts,1);

for p = 1:length(listing)
    summ = load([new_epilepsia_folder,listing(p).name]);
    summ = summ.summ;

    clinical = summ.clinical;
    %% Outcome and surgery
    all_two_year_ilae{p} = clinical.ilae{2};
    all_two_year_engel{p} = clinical.engel{2};
    all_one_year_ilae{p} = clinical.ilae{1};
    all_one_year_engel{p} = clinical.engel{1};
    all_surgery{p} = clinical.surgery;
    
end

%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    all_surgery);

%% Get outcome
switch which_outcome
    case 'ilae'
        switch which_year
            case 1
                outcome = all_one_year_ilae;
            case 2
                outcome = all_two_year_ilae;
        end
    case 'engel'
        switch which_year
            case 1
                outcome = all_one_year_engel;
            case 2
                outcome = all_two_year_engel;
        end
        
end

%% Find good and bad outcome
outcome_num = cellfun(@(x) parse_outcome(x,which_outcome),outcome);
outcome = outcome_num;

fprintf('\n%d patients had surgery and have good outcome.\n',sum(resection_or_ablation==1 & outcome==1))

end