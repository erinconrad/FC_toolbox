function intermediate_model_files

% I made this function because running the models takes a while. This
% allows me to run the models and make intermediate data files that I can
% then open for the plotting functions.

%% Parameters
pca_perc = 95;
rm_non_temporal = 1;
rm_wake = 1;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Run the lr_mt to extract features
if rm_wake == 1
    [T,features] =  lr_mt(3);
else
    error('why are you doing this?')
end

% Remove those without a response
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];

%% ROC for all features
fprintf('\nDoing main models...');
tic
just_spikes = 0;% all patients 
left = classifier_wrapper(T,features,pca_perc,1,just_spikes,rm_non_temporal,[]); % 1 means left
right = classifier_wrapper(T,features,pca_perc,2,just_spikes,rm_non_temporal,[]); % 2 means right

% Get ROC stats
[XL,YL,~,AUCL] = perfcurve(left.class,left.scores,left.pos_class);
[XR,YR,~,AUCR] = perfcurve(right.class,right.scores,right.pos_class);
left.XL = XL;
left.YL = YL;
left.AUCL = AUCL;
right.XR = XR;
right.YR = YR;
right.AUCR = AUCR;

%% ROC for spikes
just_spikes = 1; % Just spikes
lefts = classifier_wrapper(T,features,pca_perc,1,just_spikes,rm_non_temporal,[]);
rights = classifier_wrapper(T,features,pca_perc,2,just_spikes,rm_non_temporal,[]);

% Get ROC stats
[XL,YL,~,AUCL] = perfcurve(lefts.class,lefts.scores,lefts.pos_class);
[XR,YR,~,AUCR] = perfcurve(rights.class,rights.scores,rights.pos_class);

lefts.XL = XL;
lefts.YL = YL;
lefts.AUCL = AUCL;
rights.XR = XR;
rights.YR = YR;
rights.AUCR = AUCR;


%% ROC for binary spikes
just_spikes = 2; % Binary spikes
leftd = classifier_wrapper(T,features,pca_perc,1,just_spikes,rm_non_temporal,[]);
rightd = classifier_wrapper(T,features,pca_perc,2,just_spikes,rm_non_temporal,[]);

% Get ROC stats
[XL,YL,~,AUCL] = perfcurve(leftd.class,leftd.scores,leftd.pos_class);
[XR,YR,~,AUCR] = perfcurve(rightd.class,rightd.scores,rightd.pos_class);

leftd.XL = XL;
leftd.YL = YL;
leftd.AUCL = AUCL;
rightd.XR = XR;
rightd.YR = YR;
rightd.AUCR = AUCR;

%% Do subsampling analysis (this involves only spike features, only means, not SD)
fprintf('done, took %1.1f seconds',toc);
fprintf('\nStarting subsampling analysis.\n')
tic
% Run the lr_mt to extract features
[T,features,way,dur,sample,ss,durations] =  lr_mt_multitime([2 3]); 
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];

% Establish what I am varying
all_durs = unique(dur);
ndurs = length(all_durs);

all_samples = unique(sample);
nsamples = length(all_samples);

all_ss = unique(ss);
nss = length(all_ss);

% Initialize subsampling data
all_data = nan(nss,2,ndurs,nsamples); 

% Loop over nss
for iss = 1:nss
    
    curr_ss = all_ss(iss);

    all_auc_l = nan(ndurs,nsamples);
    all_auc_r = nan(ndurs,nsamples);
    
    % Loop over ndurs -  these will be different error bar points
    for id = 1:ndurs
        fprintf('\nDoing ss %d, dur %d...',iss,id);
        curr_dur = all_durs(id);
        
        
        % Loop over samples
        for is = 1:nsamples
            curr_sample = all_samples(is);

            % Get the relevant features
            relevant_features = contains(features,'_way1_') & ...
                contains(features,sprintf('_ss%d_',curr_ss)) & ...
                contains(features,sprintf('_dur%d_',curr_dur)) & ...
                contains(features,sprintf('_samp%d_',curr_sample));
            curr_features = features(relevant_features);

            % run the model
            just_spikes = 1; % Just spikes
            leftss = classifier_wrapper(T,curr_features,pca_perc,1,just_spikes,rm_non_temporal,[]);
            rightss = classifier_wrapper(T,curr_features,pca_perc,2,just_spikes,rm_non_temporal,[]);

            % Get ROC stats
            [~,~,~,AUCL] = perfcurve(leftss.class,leftss.scores,leftss.pos_class);
            [~,~,~,AUCR] = perfcurve(rightss.class,rightss.scores,rightss.pos_class);

            all_auc_l(id,is) = AUCL;
            all_auc_r(id,is) = AUCR;

        end
        fprintf('median left AUC is %1.2f\n',nanmedian(all_auc_l(id,:)))

    end

    % save data
    all_data(iss,1,:,:) = all_auc_l;
    all_data(iss,2,:,:) = all_auc_r;

end
fprintf('\nDone with subsampling analysis, took %1.1f seconds.\n',toc)


% Put everything in an out file
out.full_model.left = left;
out.full_model.right = right;
out.spike_model.left = lefts;
out.spike_model.right = rights;
out.binary_spike_model.left = leftd;
out.binary_spike_model.right = rightd;
out.subsampling.data = all_data;
out.subsampling.durations = durations;

save([plot_folder,'models.mat'],'out')

end