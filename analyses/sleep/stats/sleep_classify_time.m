function AUC = sleep_classify_time(wake_or_sleep,duration)

%{
Train and test a LR classifier using a certain total duration of
non-contiguous sleep or wake time periods. Parameters:

wake_or_sleep : 1 or 2. 1 indicates take wake periods and 2 indicates take
sleep periods

duration: a parameters indicating the number of minute-long segments to
randomly take. With shorter durations, I am using a smaller amount of time
to get spike rates to train and test the classifier at localizing the SOZ
%}

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
rate_time = out.time_out.all_spikes;
ws_time = out.time_out.all_ws;
soz = out.bin_out.all_is_soz;
pt_idx = (1:length(rate_time))';

%% Put data into friendly format
vec_rate = [];
vec_pt_idx = [];
vec_soz = [];

%% Loop over patients
for ip = 1:length(pt_idx)
    
    curr_rate_time = rate_time{ip};
    curr_ws_time = ws_time{ip,wake_or_sleep};
    curr_soz = (soz{ip})';
    
    nsoz = sum(curr_soz==1);
    nelecs = length(curr_soz);
    if randomize_soz == 1
        fake_soz_indices = randsample(nelecs,nsoz);
        curr_soz = zeros(nelecs,1);
        curr_soz(fake_soz_indices) = 1;
    end
    
    %% Get the indices of the time segments of wake or sleep
    curr_indices = find(curr_ws_time == 1); % find the indices matching the desired state
    
    % how many segments to take (either duration or all the indices)
    if isempty(duration)
        nsegments_to_take = length(curr_indices);
    else
        nsegments_to_take = min([length(curr_indices),duration]);
    end
    
    %% Get random sample of nsegments_to_take from curr_indices
    segs = randsample(curr_indices,nsegments_to_take);
    
    %% Get the rates in those
    seg_rates = curr_rate_time(:,segs);
    
    %% Take the average across times
    avg_segs_rates = nanmean(seg_rates,2);
      
    if do_norm
        % normalize the rate across electrodes (this is so that patients
        % with higher spike rates in general will get the same weight as
        % patients with lower spike rates)
        avg_segs_rates = (avg_segs_rates - nanmean(avg_segs_rates,1))./...
            nanstd(avg_segs_rates,[],1);

    end
    
    vec_rate = [vec_rate;avg_segs_rates];
    vec_pt_idx = [vec_pt_idx;repmat(pt_idx(ip),length(avg_segs_rates),1)];  
    vec_soz = [vec_soz;curr_soz];
    
end



%% Make table
T = table(vec_soz,vec_pt_idx,vec_rate);

%% Remove nan and inf rows
nan_rows = isnan(T.vec_rate) | isnan(T.vec_soz) | isnan(T.vec_pt_idx);
T(nan_rows,:) = [];


%% Divide into training and testing set
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

%% Train model with glm model
glm = fitglm(T_train,...
    'vec_soz ~ vec_rate ',...
    'Distribution','Binomial');

%% Test model on testing data
classification = predict(glm,T_test);

all_soz = T_test.vec_soz==1;
all_no_soz = T_test.vec_soz==0;

labels = cell(length(classification),1);
labels(all_soz) = {'SOZ'};
labels(all_no_soz) = {'Not SOZ'};
posclass = 'SOZ';
scores = classification;
[~,~,~,AUC,~] = perfcurve(labels,scores,posclass);

end

