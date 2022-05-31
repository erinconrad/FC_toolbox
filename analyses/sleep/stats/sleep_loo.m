function [pt_stats,X,Y,pt_specific] = sleep_loo

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
pt_stats = nan(npts,6);
pt_specific = cell(npts,3);
all_roc = cell(npts,1);
for ip = 1:npts
    fprintf('\npatient = %d of %d\n',ip,npts);
    curr_soz = all_soz{ip};
    if sum(curr_soz) == 0
        continue;
    else
        mout = updated_classifier_may2022(ip,1,[],[]);
        if ~isfield(mout,'labels'), continue; end
        for ic = 1:4
            pt_stats(ip,ic) = mout.model.Coefficients{ic+1,2}; 
        end
        
        curr_roc = [mout.X,mout.Y];
        all_roc{ip} = curr_roc;
        
        pt_stats(ip,5) = mout.AUC;
        %pt_stats(ip,6) = mout.PPV;
        %pt_stats(ip,7) = mout.NPV;
        
        pt_specific{ip,1} = mout.scores;
        pt_specific{ip,2} = mout.T;
        pt_specific{ip,3} = mout.all_soz;
        
    end
end

%% Unify roc
[X,Y] = unify_roc(all_roc);


end