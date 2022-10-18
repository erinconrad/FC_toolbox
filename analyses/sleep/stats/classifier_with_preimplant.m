function mout = classifier_with_preimplant(leave_out,wake_or_sleep,duration,just_gray,only_good_outcome,multi_stage)

%{
This is the main function to do models for the spikes and sleep paper.
Input parameters:
- leave_out: which patient to keep as testing data. This is empty if doign
bootstrap
- do_glme: do mixed effects model? 1 for paper
- wake_or_sleep: wake or sleep state. Empty if all states.
- duration of interictal data to take. Empty if full duration
- just_gray: only include gray matter electrodes?
%}

%% Parameters
do_norm = 1; % normalize spike rates within patient?
randomize_soz = 0; % negative control. If I shuffle the SOZ, I should not exceed chance AUC

%% File locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% pt file as well
pt = load([out_folder1,'pt.mat']);
pt = pt.pt;

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

% surgery with good outcome or bad outcome
surg_good = resection_or_ablation & (outcome == 1);


%% Get stuff
rate_sw = out.bin_out.all_elecs_rates_sw; %rates wake and sleep
soz = out.bin_out.all_is_soz;
pt_idx = (1:length(rate_sw))';
rate_pre_post = out.sz_out.pre_post_ictal_rates;
rate_time = out.time_out.all_spikes;
ws_time = out.time_out.all_ws;
elec_locs = out.circ_out.all_elec_locs;
names = out.circ_out.names;
npts = length(names);
ss_times = out.seeg_out.all_ss_times;

%% Get preimplant data
mri_lesional = cell(npts,1);
concordant_loc = cell(npts,1);
concordant_lat = cell(npts,1);
for ip = 1:npts
    cname = names{ip};
    found_it = 0;
    for jp = 1:length(pt)
        if strcmp(pt(jp).name,cname)
            found_it = 1;
            mri_lesional{ip} = pt(jp).clinical.pre_implant.MRI_lesional;
            concordant_loc{ip} = pt(jp).clinical.pre_implant.concordant_loc;
            concordant_lat{ip} = pt(jp).clinical.pre_implant.concordant_lat;
            break
        end
    end
    assert(found_it == 1)
end

% clean them (turn to 1s and 0s)
mri_lesional = cellfun(@(x) clean_preimplant_designations(x),mri_lesional);
concordant_loc = cellfun(@(x) clean_preimplant_designations(x),concordant_loc);
concordant_lat = cellfun(@(x) clean_preimplant_designations(x),concordant_lat);
%table(names,mri_lesional,concordant_loc,concordant_lat)

%% Put data into single vectors
vec_rate = [];
vec_rate_sleep = [];
vec_rate_wake = [];
vec_rate_post = [];
vec_rate_pre = [];
vec_pt_idx = [];
vec_soz = [];
vec_mri_lesional = [];
vec_concordant_loc = [];
vec_concordant_lat = [];

%% Which patients to do
if only_good_outcome
    pt_idx = find(surg_good);
else
    pt_idx = (1:length(rate_sw))';
end


% Loop over patients
for i = 1:length(pt_idx)
    
    ip = pt_idx(i);
    
    %% Get spike rates in different states
    curr_rate_sw = rate_sw{ip};
    curr_rate_post = rate_pre_post{ip}(:,2);
    curr_rate_pre = rate_pre_post{ip}(:,1);
    curr_soz = (soz{ip})';
    curr_rate_time = rate_time{ip};
    curr_locs = elec_locs{ip};

    
    %% remove non-gray matter?
    % Get whether electrodes are gray matter
    is_gray = strcmp(curr_locs,'other cortex') | strcmp(curr_locs,'temporal neocortical') | strcmp(curr_locs,'mesial temporal');
    
    % remove those that aren't gray matter
    if just_gray
        curr_rate_sw(~is_gray,:) = [];
        curr_rate_post(~is_gray) = [];
        curr_rate_pre(~is_gray) = [];
        curr_soz(~is_gray) = [];
        curr_rate_time(~is_gray,:) = [];
    end
    
    %% Get the indices of the time segments of wake or sleep
    if ~isempty(wake_or_sleep)
        if multi_stage
            curr_ws_time = ss_times{ip,wake_or_sleep};
        else
            curr_ws_time = ws_time{ip,wake_or_sleep};
        end
        curr_indices = find(curr_ws_time == 1); % find the indices matching the desired state
    else
        curr_indices = 1:size(curr_rate_time,2);
    end
        
    nsoz = sum(curr_soz==1);
    nelecs = length(curr_soz);
    if randomize_soz == 1
        fake_soz_indices = randsample(nelecs,nsoz);
        curr_soz = zeros(nelecs,1);
        curr_soz(fake_soz_indices) = 1;
    end
    
    % how many segments to take (either duration or all the indices)
    if isempty(duration)
        nsegments_to_take = length(curr_indices);
    else
        nsegments_to_take = min([length(curr_indices),duration]);
    end
    
    %% Get random sample of nsegments_to_take from curr_indices
    % fix for weird thing
    if length(curr_indices) == 1
        segs = curr_indices;
    else
        segs = randsample(curr_indices,nsegments_to_take);
    end
    
    %% Get the rates in those
    seg_rates = curr_rate_time(:,segs);
    
    %% Take the average across times
    avg_segs_rates = nanmean(seg_rates,2);
    if ~isempty(wake_or_sleep)
        assert(all(curr_ws_time(segs)==1))
    end
    
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
        avg_segs_rates = (avg_segs_rates - nanmean(avg_segs_rates,1))./...
            nanstd(avg_segs_rates,[],1);
    end
    
    vec_rate_sleep = [vec_rate_sleep;curr_rate_sw(:,2)];
    vec_rate_wake = [vec_rate_wake;curr_rate_sw(:,1)];
    vec_rate_pre = [vec_rate_pre;curr_rate_pre];
    vec_rate_post = [vec_rate_post;curr_rate_post];
    
    vec_pt_idx = [vec_pt_idx;repmat(ip,length(curr_rate_post),1)];

    vec_soz = [vec_soz;curr_soz];
    vec_rate = [vec_rate;avg_segs_rates];
    
    vec_mri_lesional = [vec_mri_lesional;repmat(mri_lesional(ip),length(curr_rate_post),1)];
    vec_concordant_loc = [vec_concordant_loc;repmat(concordant_loc(ip),length(curr_rate_post),1)];
    vec_concordant_lat = [vec_concordant_lat;repmat(concordant_lat(ip),length(curr_rate_post),1)];
    
end

%% Make table
T = table(vec_soz,vec_pt_idx,vec_rate_sleep,vec_rate_wake,...
    vec_rate_pre,vec_rate_post,vec_rate,...
    vec_mri_lesional,vec_concordant_loc,vec_concordant_lat);

if only_good_outcome
    assert(isequal(find(surg_good),unique(T.vec_pt_idx)))
end

%% Remove nan and inf rows
if isempty(wake_or_sleep)
    nan_rows = isnan(T.vec_rate_sleep) | isnan(T.vec_rate_wake) | isnan(T.vec_soz) ...
        | isnan(T.vec_pt_idx) | isnan(T.vec_rate_pre) | isnan(T.vec_rate_post);
else
    nan_rows = isnan(T.vec_rate)| isnan(T.vec_soz) | isnan(T.vec_pt_idx);
end
T(nan_rows,:) = [];

% shouldn't be any nans for preimplant data
assert(~any(isnan(T.vec_mri_lesional))&&~any(isnan(T.vec_concordant_loc))...
    &&~any(isnan(T.vec_concordant_lat)))

%% Divide into training and testing data
% Are we doing LOO and if so which
mout.funny_error = 0;
mout.empty_testing = 0;
mout.one_label = 0;

if isempty(leave_out)
    % if this is empty, I am doing all patients with a bootstrap sample
    % with replacement to get bootstrap CI on coefficients
    training = randsample(length(pt_idx),length(pt_idx),true);
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

    if isempty(T_test)
        mout.AUC = nan;
        mout.empty_testing = 1;
        return
    end
    % Confirm I got the right pt to leave out
    train_pts = unique(T_train.vec_pt_idx);
    test_pts = unique(T_test.vec_pt_idx);
    assert(test_pts == leave_out);
    assert(~ismember(leave_out,train_pts));
    
end

%% Train model with glm(e) model

% Designate categorical variables (patient index and preimplant variables)
T_train.vec_pt_idx = nominal(T_train.vec_pt_idx);
T_test.vec_pt_idx = nominal(T_test.vec_pt_idx);
T_train.vec_mri_lesional = nominal(T_train.vec_mri_lesional);
T_test.vec_mri_lesional = nominal(T_test.vec_mri_lesional);
T_train.vec_concordant_loc = nominal(T_train.vec_concordant_loc);
T_test.vec_concordant_loc = nominal(T_test.vec_concordant_loc);
T_train.vec_concordant_lat = nominal(T_train.vec_concordant_lat);
T_test.vec_concordant_lat = nominal(T_test.vec_concordant_lat);

if isempty(wake_or_sleep)    

    try
        %{
        glm = fitglme(T_train,...
            ['vec_soz ~ vec_rate_wake + vec_rate_sleep + '...
            'vec_rate_pre + vec_rate_post + '...
            '(1|vec_pt_idx)'],...
            'Distribution','Binomial');
        %}
        
        % Mixed effects model. Fixed effects are spike rate in each of the
        % 4 conditions as well as the three binary preimplant variables.
        % Random effect is patient.
        glm = fitglme(T_train,...
            ['vec_soz ~ vec_rate_wake + vec_rate_sleep + '...
            'vec_rate_pre + vec_rate_post + '...
            ' vec_mri_lesional + vec_concordant_loc + vec_concordant_lat +'...
            '(1|vec_pt_idx)'],...
            'Distribution','Binomial');
        %}
    catch ME

        if contains(ME.message,'NaN or Inf values are not allowed in X.')
               mout.AUC = nan;
               mout.funny_error = 1;
               return
        else
            error('what');
        end

    end
    
    
else
    
    
    

    try
        % just the model with the rate in that specific state
        glm = fitglme(T_train,...
            'vec_soz ~ vec_rate + (1|vec_pt_idx)',...
            'Distribution','Binomial');

    catch ME
        if contains(ME.message,'NaN or Inf values are not allowed in X.')
               mout.AUC = nan;
               mout.funny_error = 1;
               return
        else
            error('what');
        end
    end
    
    
end
    
%% Test model on testing data
% probability from LR model
classification = predict(glm,T_test);
scores = classification;

% True labels
all_soz = T_test.vec_soz==1;
all_no_soz = T_test.vec_soz==0;
labels = cell(length(classification),1);
labels(all_soz) = {'SOZ'};
labels(all_no_soz) = {'Not SOZ'};
posclass = 'SOZ';

% Skip if only one class in testing data
if length(unique(labels)) == 1
    mout.AUC = nan;
    mout.one_label = 1;
    return
end
[X,Y,T,AUC,~] = perfcurve(labels,scores,posclass);

if isnan(AUC) && mout.funny_error == 0 && mout.one_label == 0 && mout.empty_testing == 0
    error('why');
end

%% Also get PPV and NPV for a specific spot on curve
% Pick the threshold to be that such that the model estimates that p of the
% electrodes are the SOZ, where p is the true proportion of SOZ electrodes
%{
nsoz = sum(all_soz);

% Calculate the number of predicted SOZ for each possible threshold in T
n_pred_soz = nan(length(T),1);
for it = 1:length(T)
    n_pred_soz(it) = sum(classification >= T(it));
end
%}


if 0
   table(predicted_labels,labels,T_test.vec_rate_sleep) 
end

%% Output
mout.AUC = AUC;
mout.model = glm;
mout.X = X;
mout.Y = Y;
mout.T = T;
mout.nsoz = nsoz;
mout.nelecs = length(all_soz);
mout.scores = scores;
mout.labels = labels;
mout.all_soz = all_soz;


end