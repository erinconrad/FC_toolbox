function pt_stats = sleep_loo

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
pt_stats = nan(npts,3);
for ip = 1:npts
    fprintf('\npatient = %d of %d\n',ip,npts);
    curr_soz = all_soz{ip};
    if sum(curr_soz) == 0
        continue;
    else
        mout = updated_classifier_may2022(ip,1);
        if ~isfield(mout,'PPV'), continue; end
        pt_stats(ip,1) = mout.AUC;
        pt_stats(ip,2) = mout.PPV;
        pt_stats(ip,3) = mout.NPV;
    end
end



end