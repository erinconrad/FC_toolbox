function raster_combine_szs(pc,f)

%% Parameters
skip = 6; % show every skip electrode labels
m = 2; % montage (CAR)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/peri-ictal/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Put runs together
out = peri_ictal_grouping(pc);
name = pc.name;
fname = pc.file(f).name;
nszs = length(pc.file(f).sz);


%% Group close seizures
which_group = ones(nszs,1);
for is = 2:nszs
    sz_first_run = pc.file(f).sz(is).run(1).run_times(1);
    last_sz_last_run = pc.file(f).sz(is-1).run(end).run_times(2);
    
    if sz_first_run < last_sz_last_run
        which_group(is) = which_group(is-1);
    else
        which_group(is) = which_group(is-1) + 1;
    end
end

%% New figure per group
ngroups = length(unique(which_group));
for ig = 1:ngroups
    curr_sz = find(which_group == ig);
    
    all_net = [];
    all_ad = [];
    all_spikes = [];
    all_szs = [];
    sz_count = 0;
    
    % Loop over szs in this group
    for is = 1:length(curr_sz)
        s = curr_sz(is);
        data = out.file(f).sz(s).montage(m);
        spikes = data.spikes;
        net = data.net;
        ad = data.ad;
        labels = data.labels;
        nruns = size(spikes,2);
        sz = nruns/2;
        
        % If it's not the last one, discard the data after the sz (to
        % stitch together)
        if is < length(curr_sz)
            net = net(:,1:sz);
            ad = ad(:,1:sz);
            spikes = spikes(:,1:sz);
        end
        
        % Stitch with other szs in group
        all_net = [all_net,net];
        all_ad = [all_ad,ad];
        all_spikes = [all_spikes,spikes];
        all_szs = [all_szs;sz+count];
        count = count+nruns;
        
    end
    
    
end

end