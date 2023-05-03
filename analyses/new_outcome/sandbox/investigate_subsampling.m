function investigate_subsampling

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
edf_path = [results_folder,'edf_summ_out/'];
sleep_stage_path = [results_folder,'edf_out/'];
out_folder = [results_folder,'analysis/new_outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

%% Load out file
out = load([out_folder,'mt_out.mat']);
out = out.out;

spikes_subsample = out.spikes_subsample;

% choose a montage, sleep stage, way
iw = 1; % random
iss = 3; % sleep
im = 1; 
a = squeeze(spikes_subsample(:,im,iss,iw,:,:));

a_agreement = nan(size(a,1),1);
for i = 1:size(a,1)
    % Get one minute duration and longest duration
    short = squeeze(a(i,1,:));
    long = squeeze(a(i,5,:));

    all_agg = nan(10,1);
    for r = 1:10
        short_curr = short{r};
        long_curr = long{r};

        if ~isempty(short_curr)
            all_agg(r) = corr(short_curr,long_curr,'rows','pairwise');
        end
    end

    a_agreement(i) = nanmedian(all_agg);

end

end