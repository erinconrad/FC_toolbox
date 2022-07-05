function all_out = laterality_model(hemi_fc_lr,hemi_spikes_lr,soz_lats)

nb = 50;
prop_train = 2/3;

%% Define 3 models
models = {'fc','spikes','all'};
nmodels = length(models);

%% Remove bilateral
non_unilateral = ~strcmp(soz_lats,'right') & ~strcmp(soz_lats,'left');
soz_lats(non_unilateral) = [];
hemi_fc_lr(non_unilateral,:) = [];
hemi_spikes_lr(non_unilateral,:) = [];

%% Also remove patients with nans for spikes and fc 
any_nan = any(isnan(hemi_spikes_lr),2) & any(isnan(hemi_fc_lr),2);
soz_lats(any_nan) = [];
hemi_fc_lr(any_nan,:) = [];
hemi_spikes_lr(any_nan,:) = [];

%% Convert lats to binary: L = 1, R = 0
npts = length(soz_lats);
all_pts = 1:npts;
soz_bin = nan(npts,1);
soz_bin(strcmp(soz_lats,'left')) = 1;
soz_bin(strcmp(soz_lats,'right')) = 0;
assert(~any(isnan(soz_bin)))

%% Break into L and R
fc_l = hemi_fc_lr(:,1);
fc_r = hemi_fc_lr(:,2);
spikes_l = hemi_spikes_lr(:,1);
spikes_r = hemi_spikes_lr(:,2);

%% Make table
T = table(soz_bin,fc_l,fc_r,spikes_l,spikes_r);
all_ppv = nan(nmodels,nb);
all_npv = nan(nmodels,nb);
all_auc = nan(nmodels,nb);
all_roc = cell(nmodels,nb);

%% Loop over models

for im = 1:nmodels
    curr_model = models{im};
    more_nans = zeros(npts,1);
    currT = T;
    
    
    if im == nmodels
        coeff_estimates = nan(4,nb);
    end
    
    %% Define model formula
    switch curr_model
        case 'fc'
            formula = 'soz_bin ~ fc_l + fc_r';
            more_nans(any(isnan(hemi_fc_lr),2)) = 1;
        case 'spikes'
            formula = 'soz_bin ~ spikes_l + spikes_r';
            more_nans(any(isnan(hemi_spikes_lr),2)) = 1;
        case 'all'
            formula = 'soz_bin ~ fc_l + fc_r + spikes_l + spikes_r';
            more_nans(any(isnan(hemi_spikes_lr),2) | any(isnan(hemi_fc_lr),2)) = 1;
    end
    more_nans = logical(more_nans);
    currT(more_nans,:) = [];
    curr_npts = size(currT,1);
    curr_all_pts = 1:curr_npts;
    
        
    %% Loop over bootstraps
    ib = 1;
    while ib < nb+1
        

        %% Identify training and testing patients
        training = randsample(curr_npts,floor(curr_npts*prop_train));
        testing = all_pts(~ismember(curr_all_pts,training));
        assert(isequal(union(training,testing),curr_all_pts'))
        assert(isempty(intersect(training,testing)))

        %% Separate table
        T_train = currT(training,:);
        T_test = currT(testing,:);

        %% Train model
        glm = fitglm(T_train,formula,'Distribution','Binomial');

        %% Test model
        classification = predict(glm,T_test);

        % Redo if all testing is the same SOZ laterality
        if all(T_test.soz_bin == 1) || all(T_test.soz_bin==0)
            continue
        end
    
        
        %% Compare classification against true
        [X,Y,~,AUC,~] = perfcurve(T_test.soz_bin,classification,1);
        
        % Redo if AUC is nan
        if isnan(AUC)
            continue
        end

        %% Assign threshold of 0.5 and build confusion matrix
        ntest = size(T_test,1);
        predicted = cell(ntest,1);
        predicted(classification>=0.5) = {'left'};
        predicted(classification<0.5) = {'right'};
        true_lat = cell(ntest,1);
        true_lat(T_test.soz_bin==1) = {'left'};
        true_lat(T_test.soz_bin==0) = {'right'};
        conf = confusion_matrix(predicted,true_lat,0);
        all_ppv(im,ib) = conf.ppv;
        all_npv(im,ib) = conf.npv;
        all_auc(im,ib) = AUC;
        all_roc{im,ib} = [X,Y];
        ib = ib + 1;
        
        
        %% Also do bootstrapping to get estimates of coefficients
        if im == nmodels
            bootstrap_pts = randsample(curr_npts,curr_npts,true); % true means sample with replacement
            T_boot = currT(bootstrap_pts,:);
            glm = fitglm(T_boot,formula,'Distribution','Binomial');
            for ic = 1:4
                beta = glm.Coefficients{ic+1,1};
                coeff_estimates(ic,ib) = beta;
            end
        end
        
    end
    
    %% Get bootstrap CI and p-value for coefficients
    if im == nmodels
        coeff_mean_ci_p = nan(4,4);
        coeff_median = nan(4,1);
        for ic = 1:4
            out = bootstrap_ci_and_p(coeff_estimates(ic,:));
            coeff_mean_ci_p(ic,:) = [out.mean,out.CI_95,out.p];
            coeff_median = nanmedian(coeff_estimates(ic,:));
        end
    end
    
    
end

%% for each pair of models, generate stats for AUC difference 
model_p = nan(nmodels,nmodels);
for im = 1:nmodels
    for jm = 1:im-1
        
        % Get aucs of each
        auc_i = all_auc(im,:);
        auc_j = all_auc(jm,:);
        
        % bootstrap stats
        stats = bootstrap_ci_and_p(auc_i,auc_j);
        model_p(im,jm) = stats.p;
        
    end
end

%% Also estimate bootstap CI of own AUC
bootci = nan(nmodels,2);
for im = 1:nmodels
    auc = all_auc(im,:);
    bootci(im,:) = [prctile(auc,2.5) prctile(auc,97.5)];
end

%% Unify ROC
for im = 1:nmodels
    curr_roc = all_roc(im,:);
    [unify_x,unify_y] = unify_roc(curr_roc);
    model_info(im).x = unify_x;
    model_info(im).y = unify_y;
    model_info(im).ym = mean(unify_y);
    model_info(im).ly = prctile(unify_y,5,1);
    model_info(im).uy = prctile(unify_y,95,1);
    model_info(im).ppv = all_ppv(im,:);
    model_info(im).npv = all_npv(im,:);
    model_info(im).auc = all_auc(im,:);
end

all_out.model_info = model_info;
all_out.model_p = model_p;
all_out.bootci = bootci;
all_out.coeff_bootstrap = coeff_mean_ci_p;
all_out.coeff_median = coeff_median;

end