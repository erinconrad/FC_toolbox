function [pt_stats,X,Y,pt_specific,excluded] = sleep_loo(just_gray)

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


%% LOO analysis - get distribution of individual patient AUCs
pt_stats = nan(npts,6); % first 4 are the model coefficients, 5 is AUC, 6 doesn't exist
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
    
    if sum(curr_soz) == 0
        excluded(ip,1) = 1;
        continue;
    else
        mout = updated_classifier_may2022(ip,1,[],[],just_gray);
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
        for ic = 1:4
            pt_stats(ip,ic) = mout.model.Coefficients{ic+1,2}; 
        end
        
        
        % model ROC and other info
        curr_roc = [mout.X,mout.Y];
        all_roc{ip} = curr_roc;

        
        pt_stats(ip,5) = mout.AUC;
       % pt_stats(ip,6) = mout.PPV;
       % pt_stats(ip,7) = mout.NPV;
        
        pt_specific{ip,1} = mout.scores;
        pt_specific{ip,2} = mout.T;
        pt_specific{ip,3} = mout.all_soz;
        
    end
end

%% Unify roc (single set of Xs and then the corresponding Ys for each model)
[X,Y] = unify_roc(all_roc);


end