function all_out = laterality_loo(hemi_fc_lr,hemi_spikes_lr,soz_lats)

%% Define 3 models
models = {'fc','spikes','all'};
nmodels = length(models);

%% Remove bilateral
removed.noriginal = length(soz_lats);
non_unilateral = ~strcmp(soz_lats,'right') & ~strcmp(soz_lats,'left');
soz_lats(non_unilateral) = [];
hemi_fc_lr(non_unilateral,:) = [];
hemi_spikes_lr(non_unilateral,:) = [];
removed.nremoved_non_unilateral = sum(non_unilateral);


%% Also remove patients with nans for spikes and fc 
removed.nremoved_missing_spikes = sum(any(isnan(hemi_spikes_lr),2));
removed.nremoved_missing_fc = sum(any(isnan(hemi_fc_lr),2));
any_nan = any(isnan(hemi_spikes_lr),2) & any(isnan(hemi_fc_lr),2);
soz_lats(any_nan) = [];
hemi_fc_lr(any_nan,:) = [];
hemi_spikes_lr(any_nan,:) = [];


%% Convert lats to binary: L = 1, R = 0
npts = length(soz_lats);
soz_bin = nan(npts,1);
soz_bin(strcmp(soz_lats,'left')) = 1;
soz_bin(strcmp(soz_lats,'right')) = 0;
assert(~any(isnan(soz_bin)))

%% Break into L and R
fc_l = hemi_fc_lr(:,1);
fc_r = hemi_fc_lr(:,2);
spikes_l = hemi_spikes_lr(:,1);
spikes_r = hemi_spikes_lr(:,2);

%% Make patient identifier
new_npts = length(fc_l);
pt_id = (1:new_npts)';

%% Make table
T = table(pt_id,soz_bin,fc_l,fc_r,spikes_l,spikes_r);


%% Loop over models
conf_table_labels = {'left-correct','left-incorrect';'right-incorrect','right-correct'};
all_conf_table = zeros(2,2,nmodels);
all_lpv = zeros(nmodels,1);
all_rpv = zeros(nmodels,1);
all_accuracy = zeros(nmodels,1);
nremoved = nan(nmodels,1);
nfinal = nan(nmodels,1);

for im = 1:nmodels
    curr_model = models{im};
    more_nans = zeros(npts,1);
    currT = T;
    conf_table = zeros(2,2);
   
    
    %% Define model formula
    switch curr_model
        case 'fc'
            formula = 'soz_bin ~ fc_l + fc_r';
            more_nans(any(isnan(hemi_fc_lr),2)) = 1;
        case 'spikes'
            formula = 'soz_bin ~ spikes_l + spikes_r';
            more_nans(any(isnan(hemi_spikes_lr),2) | any(isnan(hemi_fc_lr),2)) = 1; % also remove absent FC to have same number
        case 'all'
            formula = 'soz_bin ~ fc_l + fc_r + spikes_l + spikes_r';
            more_nans(any(isnan(hemi_spikes_lr),2) | any(isnan(hemi_fc_lr),2)) = 1;
    end
    more_nans = logical(more_nans);
    currT(more_nans,:) = [];
    curr_npts = size(currT,1);
    curr_all_pts = 1:curr_npts;
    nremoved(im) = npts-curr_npts;
    
        
    %% Loop over patients
    for ip = 1:curr_npts
        

        %% Identify training and testing patients
        testing = ip;
        training = curr_all_pts(~ismember(curr_all_pts,testing));
        assert(isequal(union(training,testing),curr_all_pts))
        assert(isempty(intersect(training,testing)))

        %% Separate table
        T_train = currT(training,:);
        T_test = currT(testing,:);

        %% Train model
        glm = fitglm(T_train,formula,'Distribution','Binomial');

        %% Test model
        classification = predict(glm,T_test);

        %% Assign threshold of 0.5 and build confusion matrix
        if classification >=0.5
            predicted_soz = 'left';
        elseif classification < 0.5
            predicted_soz = 'right';
        end
        
        if T_test.soz_bin==1
            true_soz = 'left';
        elseif T_test.soz_bin==0
            true_soz = 'right';
        end
        
        % Add patient to correct cell of confusion matrix
        if strcmp(predicted_soz,'left') && strcmp(true_soz,'left')
            conf_table(1,1) = conf_table(1,1)+1;
        elseif strcmp(predicted_soz,'right') && strcmp(true_soz,'left')
            conf_table(1,2) = conf_table(1,2)+1;
        elseif strcmp(predicted_soz,'left') && strcmp(true_soz,'right')
            conf_table(2,1) = conf_table(2,1)+1;
        elseif strcmp(predicted_soz,'right') && strcmp(true_soz,'right')
            conf_table(2,2) = conf_table(2,2)+1;
        else
            error('why');
        end
        

    end
    
    lpv = conf_table(1,1)/(conf_table(1,1)+conf_table(2,1));
    rpv = conf_table(2,2)/(conf_table(2,2)+conf_table(1,2));
    accuracy = (conf_table(1,1)+conf_table(2,2))/...
        (conf_table(1,1)+conf_table(2,2)+conf_table(1,2)+conf_table(2,1));
    all_conf_table(:,:,im) =  conf_table;
    all_lpv(im) = lpv;
    all_rpv(im) = rpv;
    all_accuracy(im) = accuracy;
    nfinal(im) = sum(conf_table(1,1)+conf_table(2,2)+conf_table(1,2)+conf_table(2,1));
    
end


removed.nremoved = nremoved;
removed.nfinal = nfinal;

all_out.removed = removed;
all_out.all_conf_table = all_conf_table;
all_out.conf_table_labels = conf_table_labels;
all_out.all_lpv = all_lpv;
all_out.all_rpv = all_rpv;
all_out.all_accuracy = all_accuracy;


end