function [pt_stats,X,Y,pt_specific,excluded] = sleep_loo(just_gray,only_good_outcome)

%% Locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];
time_roc_folder = [out_folder,'time_roc/'];
if ~exist(time_roc_folder,'dir')
    mkdir(time_roc_folder)
end

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;
npts = length(out.circ_out.names);
all_soz = out.bin_out.all_is_soz;

%% Get some outcome stuff
two_year_engel = out.circ_out.all_two_year_engel;
two_year_ilae = out.circ_out.all_two_year_ilae;
surgery = out.circ_out.all_surgery;

% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);

% Parse outcome
outcome = cellfun(@(x) parse_outcome(x,'engel'),two_year_engel);
%outcome = cellfun(@(x) parse_outcome(x,'ilae'),two_year_ilae);

% surgery with good outcome or bad outcome
surg_good = resection_or_ablation & (outcome == 1);
surg_bad = resection_or_ablation & (outcome == 0);



%% LOO analysis - get distribution of individual patient AUCs
pt_stats = nan(npts,8); % first 7 are the model coefficients, 8 is AUC
pt_specific = cell(npts,3); % model scores, then potential thresholds, then SOZs for each patient
all_roc = cell(npts,1); % individual patient ROCs
excluded = zeros(npts,4); % empty soz, empty testing, funny error, one label

%% Notes on patient exclusions
% patients 2, 42, 71: no SOZ
% patient 7 (HUP108): no sleep times identified (looking at raster,almost
% all bad times)
% patient 17 (HUP122): no postictal times, no seizures in seizure time file despite
% SOZ being listed
% patient 49: wacky error, model fails to converge

for ip = 1:npts
    %fprintf('\npatient = %d of %d\n',ip,npts);
    curr_soz = all_soz{ip};

    % Don't do if bad outcome
    if surg_good(ip) == 0 && only_good_outcome == 1
        continue
    end
    
    if sum(curr_soz) == 0
        excluded(ip,1) = 1;
        continue;
    else
        mout = classifier_with_preimplant(ip,[],[],just_gray,only_good_outcome);
        %mout = updated_classifier_may2022(ip,1,[],[],just_gray,pre_implant);
        % first argument indicates which patient to be held out as testing
        % data.
        
        if mout.empty_testing == 1
            excluded(ip,2) = 1;
        end
        
        if mout.funny_error == 1
            excluded(ip,3) = 1;
        end
        
        if mout.one_label == 1
            excluded(ip,4) = 1;
        end
        
        if ~isfield(mout,'labels'), continue; end
        
        % Get model coefficients
        for ic = 1:7
            pt_stats(ip,ic) = mout.model.Coefficients{ic+1,2}; 
        end
        
        
        % model ROC and other info
        curr_roc = [mout.X,mout.Y];
        all_roc{ip} = curr_roc;

        
        pt_stats(ip,8) = mout.AUC;
        
        pt_specific{ip,1} = mout.scores;
        pt_specific{ip,2} = mout.T;
        pt_specific{ip,3} = mout.all_soz;
        
    end
end

%% Unify roc (single set of Xs and then the corresponding Ys for each model)
[X,Y] = unify_roc(all_roc);


end