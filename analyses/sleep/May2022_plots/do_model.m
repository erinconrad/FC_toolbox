function do_model
%{
1) It looks like sleep and post-ictal spike rates independently help localize,
but awake and pre-ictal do not. (CONFIRM)
2) You only need about 30 minutes of sleep to get the most bang for your buck (the other analysis)
3) Patients are heterogeneous in how well spikes localize, with TLE
patients doing much better
4) I need to think of some way to conceptualize actual thresholds. Perhaps
I should set thresholds such that same number of SOZ as median number. How
many false positives and false negatives?? If x is the proportion of
electrodes that I want to be designated SOZ electrodes, that means I want
(TP + FP)/(TP+TP+FN+FP) =  X
I also have the constraint of the ROC curve
%}


%% Parameters
nb = 10;
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

%% Get proportion of electrodes that are soz
nsoz = cellfun(@(x) sum(x==1),all_soz);
nelecs = cellfun(@length,all_soz);
median_proportion = median(nsoz./nelecs);

%% Get some clinical stuff
stereo = logical(out.circ_out.stereo);
loc = out.circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');

%% Do main model both to get estimates of coefficients and AUCs
all_auc = nan(nb,1);
all_coeffs = nan(ncoeffs,nb);
all_roc = cell(nb,1);
for ib = 1:nb
     fprintf('\nib = %d of %d\n',ib,nb)
     mout = updated_classifier_may2022([],1);
     all_auc(ib) = mout.AUC;
     for ic = 1:ncoeffs
            % get the estimate of the model coefficient
            all_coeffs(ic,ib) = mout.model.Coefficients{ic+1,2}; 
         
     end
     all_roc{ib} = [mout.X,mout.Y];
    
end
[X,Y] = unify_roc(all_roc);

%% Calculate thresholds
% Find the point along the ROC curve that would return about 5% positivity
% rate
closest_point = find_point_roc_proportion(X,mean(Y,1),median_proportion);


% Get confusion matrices
cout = confusion_matrix(predicted,actual,do_plot);

%% Determine bootstrap stats for each coefficient
coeff_stats = nan(ncoeffs,4); % mean, lower CI, higher CI, p
for ic = 1:ncoeffs
    tout = bootstrap_ci_and_p(squeeze(all_coeffs(ic,:)));
    coeff_stats(ic,:) = [tout.mean,tout.CI_95,tout.p];
    
end


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

% Stereo doesn't matter
if 0
    figure
    plot(1+randn(sum(stereo),1)*0.05,pt_auc(stereo),'o','linewidth',2)
    hold on
    plot(2+randn(sum(~stereo),1)*0.05,pt_auc(~stereo),'o','linewidth',2)

end

% TLE vs eTLE DOES matter
if 0
    figure
    plot(1+randn(sum(temporal),1)*0.05,pt_auc(temporal),'o','linewidth',2)
    hold on
    plot(2+randn(sum(extra),1)*0.05,pt_auc(extra),'o','linewidth',2)

end





end