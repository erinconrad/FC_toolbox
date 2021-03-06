function predict_spike_spread

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
int_folder = [results_folder,'analysis/intermediate/'];
out_folder = [results_folder,'spike_spread/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

all_coaw = [];
all_fcw = [];
all_distw = [];
all_pt_idx = [];

%% Loop over patients
for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    coa = summ.coa;
    avg_fc = summ.avg_fc;
    labels = summ.labels;
    name = summ.name;
    locs = summ.locs;
    
    %% Expand coa and time average
    coa = wrap_or_unwrap_adjacency_fc_toolbox(coa);
    coa = nanmean(coa,3);
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    labels = labels(~ekg);
    avg_fc = avg_fc(~ekg,~ekg);
    coa = coa(~ekg,~ekg);
    locs = locs(~ekg,:);
    
    %% Build distance network
    dist_net = make_dist_network(locs);
    
    %% Normalize coa and FC by channel
    coa = (coa - nanmean(coa,2))./nanstd(coa,[],2);
    avg_fc = (avg_fc-nanmean(avg_fc,2))./nanstd(avg_fc,[],2);
    dist_net = (dist_net - nanmean(dist_net,2))./nanstd(dist_net,[],2);
    
    if 0
        figure
        nexttile
        imagesc(avg_fc)
        
        nexttile
        imagesc(coa)
        
        nexttile
        imagesc(dist_net)
    end
    
    coaw = wrap_or_unwrap_adjacency_fc_toolbox(coa);
    fcw = wrap_or_unwrap_adjacency_fc_toolbox(avg_fc);
    distw = wrap_or_unwrap_adjacency_fc_toolbox(dist_net);
    pt_idx = p*ones(length(coaw),1);
    
    all_coaw = [all_coaw;coaw];
    all_fcw = [all_fcw;fcw];
    all_distw = [all_distw;distw];
    all_pt_idx = [all_pt_idx;pt_idx];
    
    %T = table(coaw,fcw,distw);
    %mdl = fitlm(T,'coaw~fcw+distw');
    
end


T = table(all_coaw,all_fcw,all_distw,all_pt_idx);
out.T = T;
save([out_folder,'out.mat'],'out');

end
