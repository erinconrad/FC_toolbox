function do_model

%% Parameters
nb = 50;
ncoeffs = 4;
myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];

%% Seed rng (for splitting testing and training data)
% If I don't seed this the AUC usually bounces around 0.73-0.82
rng(0)

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
pt_auc = nan(npts,1);
for ip = 1:npts
    fprintf('\npatient = %d of %d\n',ip,npts);
    curr_soz = all_soz{ip};
    if sum(curr_soz) == 0
        pt_auc(ip) = nan;
    else
        mout = updated_classifier_may2022(ip,1);
        pt_auc(ip) = mout.AUC;
    end
end

%% Get some clinical stuff

%% Do main model both to get estimates of coefficients and AUCs
all_auc = nan(nb,1);
all_coeffs = nan(ncoeffs,nb);
for ib = 1:nb
     fprintf('\nib = %d of %d\n',ib,nb)
     mout = updated_classifier_may2022([],1);
     all_auc(ib) = mout.AUC;
     for ic = 1:ncoeffs
            % get the estimate of the model coefficient
            all_coeffs(ic,ib) = mout.model.Coefficients{ic+1,2}; 
         
     end
    
end

%% Determine bootstrap stats for each coefficient
coeff_stats = nan(ncoeffs,4); % mean, lower CI, higher CI, p
for ic = 1:ncoeffs
    tout = bootstrap_ci_and_p(squeeze(all_coeffs(ic,:)));
    coeff_stats(ic,:) = [tout.mean,tout.CI_95,tout.p];
    
end



end