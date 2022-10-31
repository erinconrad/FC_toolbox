function summ = sw_ad_erin_designations

%{
This function takes Erin's sleep/wake designations from the
erin_designations folder and decides whether it is a sleep state or an
awake state, and grabs the corresponding alpha delta ratio, as well as all
alpha delta ratios (for the purpose of normalization).
%}

%% Parameters
sleep_names = {'1','2','3'}; % these labels mean sleep
wake_names = {'a'}; % this label means awake
uncertain_names = {'u'};

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
erin_des_folder = [results_folder,'analysis/sleep/erin_designations/'];
spikes_folder = [results_folder,'all_out/'];
%data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Loop through all patients in the erin_designations folder
listing = dir([erin_des_folder,'*.mat']);


for l = 1:length(listing)
    
    fname = listing(l).name;
    name = strrep(fname,'.mat','');
    
    % REMOVE THIS ONCE PT DONE!!!
    %if strcmp(name,'HUP213'), continue; end
    
    % Load the file
    out = load([erin_des_folder,name,'.mat']);
    if isfield(out,'out')
        out = out.out;
    elseif isfield(out,'nout')
        out = out.nout;
    end
    
    % Skip if there's no corresponding spike (and more relevantly alpha delta) file
    if ~exist([spikes_folder,name,'_pc.mat']), continue; end
    
    % Load the spikes/ad file
    pc = load([spikes_folder,name,'_pc.mat']);
    pc = pc.pc;

    % Initialize the sleep state designations and an array of alpha delta
    % ratios corresponding to this
    sleep_state = {};
    ad = [];
    
    % Loop over the blocks I have designations for
    for f = 1:length(pc.file)
        blocks = out.file(f).blocks;
        for ib = 1:length(blocks)
            
            designation = out.file(f).block(blocks(ib)).erin;
            
            % remove non-intracranial
            labels = pc.file(f).run(blocks(ib)).data.clean_labels;
            ekg = find_non_intracranial(labels);
            
            % Get alpha delta ratio for that block
            alpha_delta = (pc.file(f).run(blocks(ib)).data.montage(2).ad);
            alpha_delta = nanmean(alpha_delta(~ekg)); % ignore ekg and scalp channels
            
            % fill these up
            sleep_state = [sleep_state;designation];
            ad = [ad;alpha_delta];
            
        end
        
    end
    
    all_ad = [];
    % Also loop over every block to get details on ad
    for f = 1:length(pc.file)
        for ib = 1:length(pc.file(f).run)
            alpha_delta = pc.file(f).run(ib).data.montage(2).ad;
            
            % remove non-intracranial
            labels = pc.file(f).run(ib).data.clean_labels;
            ekg = find_non_intracranial(labels);
            alpha_delta = (alpha_delta(~ekg));
            
            all_ad = [all_ad;nanmean(alpha_delta)];
        end
    end
    
    % find those matching my designated wake and sleep labels
    sleep_idx = ismember(sleep_state,sleep_names);
    wake_idx = ismember(sleep_state,wake_names);
    uncertain_idx = ismember(sleep_state,uncertain_names);
    
    % put the corresponding ad ratios into the right bins
    summ(l).sw.sleep = ad(sleep_idx);
    summ(l).sw.wake = ad(wake_idx);
    summ(l).sw.uncertain = ad(uncertain_idx);
    summ(l).name = name;
    summ(l).ad = all_ad;
    summ(l).designation = designation;
    
    if 0
        figure
        plot(1+randn(sum(sleep_idx),1)*0.05,ad(sleep_idx),'o')
        hold on
        plot(2+randn(sum(wake_idx),1)*0.05,ad(wake_idx),'o')
    end
    
end




end