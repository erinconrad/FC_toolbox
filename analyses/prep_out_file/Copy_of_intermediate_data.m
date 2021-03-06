function intermediate_data

%% Parameters
m = 2; % do not change
net_m = 2;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/intermediate/'];
spikes_folder = [results_folder,'all_out/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(locations.bct));

validation_file = [scripts_folder,'spike_detector/Manual validation.xlsx'];

% Pt struct
data_folder = [locations.main_folder,'data/'];
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Get the indices of the patients with good spikes
T = readtable(validation_file);
good_pts = T.Var13;
good_pts = good_pts(~isnan(good_pts));
npts = length(good_pts);

%% Also grab SOZ info
szT = readtable(validation_file,'Sheet','SOZ');

count = 0;
% Loop over patients
for l = 1:npts
    j = good_pts(l);
    name = pt(j).name;
    rid = pt(j).rid;
    
    %% Load the spike file
    fname = [spikes_folder,name,'_pc.mat'];
    if ~exist(fname,'file')
        fprintf('\nCannot find spike file for %s, skipping...\n',name);
        continue
    end
    
    %% Get block dur
    block_dur = diff(pt(j).ieeg.file(1).block_times(1,:));
    
     %% Get basic info from the patient
    % load the spike file
    pc = load(fname);
    pc = pc.pc;
    
    % Skip the patient if it's incomplete
    if length(pc.file) < length(pt(j).ieeg.file) || ...
            length(pc.file(end).run) < size(pt(j).ieeg.file(end).run_times,1)
        fprintf('\n%s incomplete, skipping\n',name);
        continue
    end
        
    % add count
    count = count + 1;
    
    % reconcile files (deal with changes in electrode names)
    out = net_over_time(pc,pt,j);
    out = reconcile_files(out);
    
    
    % Get the correct row of the SOZ table
    szr = nan;
    for k = 1:size(szT,1)
        if strcmp(szT.name{k},name)
            szr = k;
            break
        end
    end
    
    % parse the soz names
    all_labels = decompose_labels(out.all_labels,name);
    [soz_labels,soz_chs] = parse_soz(szT.SOZElectrode{k},all_labels,name);
    summ.soz.labels = soz_labels;
    summ.soz.chs = soz_chs;
    summ.soz.loc = szT.region{k};
    summ.soz.lat = szT.lateralization{k};
    
    
    % Get the spikes and the labels
    times = out.times;
    spikes = out.montage(m).spikes;
    labels = out.montage(m).labels;
    coa = out.montage(m).coa;
    rl = out.montage(m).rl;
    %coi_ch = out.montage(m).coi_ch;
    coi_global = out.montage(m).coi_global;
    file_times = out.run_center;
    file_index = out.file_index;
    ad = out.montage(m).ad;
    n_rm_ictal = out.montage(m).n_rm_ictal;
    sz_times = out.sz_times;
    seq_info = out.montage(m).seq_info;
    leader = out.montage(m).leader_montage;
    mod_midnight = out.mod_midnight;
    
    %% Bipolar networks
    % car montage for networks
    fc_bi = out.montage(1).net;
    fc_bi = wrap_or_unwrap_adjacency_fc_toolbox(fc_bi);
    
    % Get global efficiency and avg node strength over time
    %ge = ge_over_blocks(fc); % this takes a long time.
    ns_bi = ns_over_blocks(fc_bi);
    
    % Get avg fc
    avg_fc_bi = nanmean(fc_bi,3);
    
    %% CAR networks
    % car montage for networks
    fc_car = out.montage(2).net;
    fc_car = wrap_or_unwrap_adjacency_fc_toolbox(fc_car);
    
    % Get global efficiency and avg node strength over time
    %ge = ge_over_blocks(fc); % this takes a long time.
    ns_car = ns_over_blocks(fc_car);
    
    % Get avg fc
    avg_fc_car = nanmean(fc_car,3);
    
    %% make main metrics
    if net_m == 1
        ns = ns_bi;
        avg_fc = avg_fc_bi;
    elseif net_m == 2
        ns = ns_car;
        avg_fc = avg_fc_car;
    end
        
    
    % Clean the labels
    clean_labels = decompose_labels(labels,name);   
    
    % Get number of electrode localizations
    ne = length(pt(j).elecs);
    
    %% Find the anatomy corresponding to the spike labelsomy
    % Initialize a new spike_anatomy cell array corresponding to spike
    % labels
    spike_anatomy = cell(length(labels),1);
    spike_locs = nan(length(labels),3);
    already_filled = zeros(length(labels),1);
    
    for e = 1:ne
        % Get loc/anatomy names and labels
        ana_name = decompose_labels(pt(j).elecs(e).elec_names,name);
        locs = pt(j).elecs(e).locs;
        anatomy = pt(j).elecs(e).anatomy;

        % Indices of the loc/anatomy names that match the spike labels
        [lia,locb] = ismember(clean_labels,ana_name);
        % sanity check
        if ~isequal(clean_labels(lia~=0 & already_filled == 0),ana_name(locb(lia~=0 & already_filled == 0))), error('oh no'); end

        
        % Fill up spike anatomy and locs with the anatomy
        if ~strcmp(class(anatomy),'double')
            %spike_anatomy(lia~=0 & already_filled == 0) = anatomy(locb(lia~=0 & already_filled == 0));
            spike_anatomy(lia~=0) = anatomy(locb(lia~=0));
        end
        %spike_locs(lia~=0 & already_filled == 0,:) = locs(locb(lia~=0 & already_filled == 0),:);
        spike_locs(lia~=0,:) = locs(locb(lia~=0),:);
        
        % set already filled
        %already_filled(lia~=0 & already_filled == 0) = 1;
        already_filled(lia~=0) = 1;
    end
    
    % Skip this patient if most anatomy descriptions empty
    perc_empty = sum(cellfun(@(x) isempty(x),spike_anatomy))/length(spike_anatomy);
    if perc_empty > 0.5
        fprintf('\n%s has too many empty anatomical locations\n');
        loc = cell(length(spike_anatomy),1);
        lat = cell(length(spike_anatomy),1);
        bad_anatomy_flag = 1;
    else
        % Get anatomical groupings
        [loc,lat] = cluster_anatomical_location(spike_anatomy);
        bad_anatomy_flag = 0;
    end
    
    %% Put it all in the intermediate struct
    %{
    summ(count).name = name;
    summ(count).times = times;
    summ(count).spikes = spikes;
    %summ(count).coi_global = coi_global;
    summ(count).rl = rl;
    summ(count).coa = coa;
    summ(count).labels = clean_labels;
    summ(count).locs = spike_locs;
    summ(count).anatomy = spike_anatomy;
    summ(count).bad_anatomy_flag = bad_anatomy_flag;
    summ(count).ana_loc = loc;
    summ(count).ana_lat = lat;
    summ(count).file_times = file_times;
    summ(count).file_index = file_index;
    summ(count).ad = ad;
    summ(count).block_dur = block_dur;
    %}
    
    summ.name = name;
    summ.times = times;
    summ.spikes = spikes;
    summ.coi_global = coi_global;
    summ.rl = rl;
    %summ.coa = coa;
    summ.labels = clean_labels;
    summ.locs = spike_locs;
    summ.anatomy = spike_anatomy;
    summ.bad_anatomy_flag = bad_anatomy_flag;
    summ.ana_loc = loc;
    summ.ana_lat = lat;
    summ.file_times = file_times;
    summ.file_index = file_index;
    summ.ad = ad;
    summ.block_dur = block_dur;
    summ.ns = ns;
    summ.n_rm_ictal = n_rm_ictal;
    summ.sz_times = sz_times;
    summ.rid = rid;
    summ.avg_fc = avg_fc;
    summ.seq_info = seq_info;
    summ.clinical = pt(j).clinical;
    summ.leader = leader;
    summ.mod_midnight = mod_midnight;
    summ.ns_car = ns_car;
    summ.avg_fc_car = avg_fc_car;
    summ.ns_bi = ns_bi;
    summ.avg_fc_bi = avg_fc_bi;
    %summ.ge = ge;

    %% Save it all
    save([out_folder,name,'.mat'],'summ');

end




