function mout = updated_classifier_may2022(leave_out,do_glme)

%% Parameters
prop_train = 2/3;
do_norm = 1;
randomize_soz = 0; % negative control. If I shuffle the SOZ, I should not exceed chance AUC


locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];


%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Get stuff
rate_sw = out.bin_out.all_elecs_rates_sw; %rates wake and sleep
soz = out.bin_out.all_is_soz;
pt_idx = (1:length(rate_sw))';
rate_pre_post = out.sz_out.pre_post_ictal_rates;

%% Put data into friendly format
vec_rate_sleep = [];
vec_rate_wake = [];
vec_rate_post = [];
vec_rate_pre = [];
vec_pt_idx = [];
vec_soz = [];
for ip = 1:length(pt_idx)
    
    curr_rate_sw = rate_sw{ip};
    curr_rate_post = rate_pre_post{ip}(:,2);
    curr_rate_pre = rate_pre_post{ip}(:,1);
    
    if do_norm
        % normalize the rate across electrodes (this is so that patients
        % with higher spike rates in general will get the same weight as
        % patients with lower spike rates)
        curr_rate_sw = (curr_rate_sw - nanmean(curr_rate_sw,1))./...
            nanstd(curr_rate_sw,[],1);        
        curr_rate_post = (curr_rate_post - nanmean(curr_rate_post))./...
            nanstd(curr_rate_post);
        curr_rate_pre = (curr_rate_pre - nanmean(curr_rate_pre))./...
            nanstd(curr_rate_pre);
    end
    
    vec_rate_sleep = [vec_rate_sleep;curr_rate_sw(:,2)];
    vec_rate_wake = [vec_rate_wake;curr_rate_sw(:,1)];
    vec_rate_pre = [vec_rate_pre;curr_rate_pre];
    vec_rate_post = [vec_rate_post;curr_rate_post];
    
    vec_pt_idx = [vec_pt_idx;repmat(pt_idx(ip),length(rate_sw{ip}(:,2)),1)];
    
    vec_soz = [vec_soz;soz{ip}'];
    
end

%% Make table
T = table(vec_soz,vec_pt_idx,vec_rate_sleep,vec_rate_wake,vec_rate_pre,vec_rate_post);

%% Remove nan and inf rows
nan_rows = isnan(T.vec_rate_sleep) | isnan(T.vec_rate_wake) | isnan(T.vec_soz) ...
    | isnan(T.vec_pt_idx) | isnan(T.vec_rate_pre) | isnan(T.vec_rate_post);
T(nan_rows,:) = [];

%% Divide into training and testing data
% Are we doing LOO and if so which
if isempty(leave_out)
    training = randsample(length(pt_idx),floor(length(pt_idx)*prop_train));
    training_idx = ismember(T.vec_pt_idx,training);
    testing_idx = ~ismember(T.vec_pt_idx,training);

    T_train = T(training_idx,:);
    T_test = T(testing_idx,:);

    % Confirm that I separated by patients
    train_pts = unique(T_train.vec_pt_idx);
    test_pts = unique(T_test.vec_pt_idx);
    assert(isempty(intersect(train_pts,test_pts)))
    assert(sum(isnan(table2array(T)),'all')==0)
else
    testing = leave_out;
    training_idx = ~ismember(T.vec_pt_idx,testing);
    testing_idx = ismember(T.vec_pt_idx,testing);
    
    T_train = T(training_idx,:);
    T_test = T(testing_idx,:);

    
    % Confirm I got the right pt to leave out
    train_pts = unique(T_train.vec_pt_idx);
    test_pts = unique(T_test.vec_pt_idx);
    assert(test_pts == leave_out);
    assert(~ismember(leave_out,train_pts));
    
end

%% Train model with glm model
if do_glme
    T_train.vec_pt_idx = nominal(T_train.vec_pt_idx);
    glm = fitglme(T_train,...
        'vec_soz ~ vec_rate_wake + vec_rate_sleep + vec_rate_pre + vec_rate_post + (1|vec_pt_idx)',...
        'Distribution','Binomial');
else
    glm = fitglm(T_train,...
        'vec_soz ~ vec_rate_wake + vec_rate_sleep + vec_rate_pre + vec_rate_post',...
        'Distribution','Binomial');
end

%% Test model on testing data
if do_glme
    mout.glme = glm;
else
    classification = predict(glm,T_test);

    all_soz = T_test.vec_soz==1;
    all_no_soz = T_test.vec_soz==0;

    labels = cell(length(classification),1);
    labels(all_soz) = {'SOZ'};
    labels(all_no_soz) = {'Not SOZ'};
    posclass = 'SOZ';
    scores = classification;
    [~,~,~,AUC,~] = perfcurve(labels,scores,posclass);

    mout.AUC = AUC;
end

end