function prep_data_mt

%% Parameters
durations = [1 5 10 20 30];
samples_per_duration = 20;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
edf_path = [results_folder,'edf_summ_out_epilepsy_laterality/'];
sleep_stage_path = [results_folder,'edf_out/'];
out_folder = [results_folder,'analysis/new_outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%{
%% Load the current output file
data = load([out_folder,'main_out.mat']);
data = data.out;
all_names = data.all_names;
npts = length(all_names);
%}

%% Load pt folder to get names
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
npts = length(pt);
all_names = cell(npts,1);
for ip = 1:npts
    all_names{ip} = pt(ip).name;
end

% Clinical stuff
all_surgery = cell(npts,1);
all_resec_lat = cell(npts,1);
all_resec_loc = cell(npts,1);
all_ablate_lat = cell(npts,1);
all_ablate_loc = cell(npts,1);
all_engel = cell(npts,2);
all_ilae = cell(npts,2);
all_soz_loc = cell(npts,1);
all_soz_lat = cell(npts,1);

% Goal is to build in same structure
% Structure is which montage, then all vs wake vs sleep
all_coh = cell(npts,3,3);
all_bp = cell(npts,3,3);
all_spikes = cell(npts,3,3);
all_rl = cell(npts,3,3);
all_pearson = cell(npts,3,3);
all_pearson_squared = cell(npts,3,3);
all_plv = cell(npts,3,3);
all_labels = cell(npts,3);
all_se = cell(npts,3,3);
all_re = cell(npts,3,3);
all_xcor = cell(npts,3,3);
all_lags = cell(npts,3,3);
all_ll = cell(npts,3,3);
all_rel_bp = cell(npts,3,3);
all_coh_iqr = cell(npts,3,3);
all_bp_iqr = cell(npts,3,3);
all_spikes_iqr = cell(npts,3,3);
all_rl_iqr = cell(npts,3,3);
all_pearson_iqr = cell(npts,3,3);
all_pearson_squared_iqr = cell(npts,3,3);
all_plv_iqr = cell(npts,3,3);
all_se_iqr = cell(npts,3,3);
all_re_iqr = cell(npts,3,3);
all_xcor_iqr = cell(npts,3,3);
all_lags_iqr = cell(npts,3,3);
all_rel_bp_iqr = cell(npts,3,3);
all_ll_iqr = cell(npts,3,3);
all_n_wake_sleep = nan(npts,2);
all_disconnected = nan(npts,72);
all_n_wake_sleep_connected = nan(npts,2);
all_atropos = cell(npts,1);
all_dkt = cell(npts,1);

% Initialize subsampling stuff
spikes_subsample = cell(npts,3,3,2,length(durations),samples_per_duration);


%% Loop over patients
for p = 1:npts

    name = all_names{p};


    % Load the edf summary file
    if exist([edf_path,name,'/summ.mat'],'file') == 0, continue; end
    info = load([edf_path,name,'/summ.mat']);
    info = info.out;

    % Clinical stuff
    all_surgery{p} = info.surgery;
    all_resec_lat{p} = info.resec_lat;
    all_resec_loc{p} = info.resec_loc;
    all_ablate_lat{p} = info.ablate_lat;
    all_ablate_loc{p} = info.ablate_loc;
    if isempty(info.engel)
         info.engel = {nan nan};
    end
    if isempty(info.ilae)
        info.ilae = {nan nan};
    end
    all_engel(p,:) = info.engel;
    all_ilae(p,:) = info.ilae;
    all_soz_loc{p} = info.soz_loc;
    all_soz_lat{p} = info.soz_lat;

    % get the time varying measurements
    spike_counts = info.all_spike_counts;
    rl = info.all_rl;
    bp = info.all_bp;
    rel_bp = info.all_rel_bp;
    coh = info.all_coh;
    re = info.all_re;
    pc = info.all_pc;
    pc_squared = info.all_pc_squared;
    xcor = info.all_xcor;
    lags = info.all_lags;
    se = info.all_se;
    times = info.all_times;
    labels = info.montage_labels;
    ll = info.all_ll;
    
    fs = info.fs;
    plv = info.all_plv;

    rl  = rl/fs;

    % fill labels
    all_labels{p,1} = labels{1,1};
    all_labels{p,2} = labels{2,1};
    all_labels{p,3} = labels{3,1};
    clean_labels = cellfun(@(x) strrep(x,'-Ref',''),labels{1,1},'UniformOutput',false);

    % get atlas parcellations
    if 0
        if ~isempty(pt(p).atropos)
            aatlas = pt(p).atropos.label;
            alabels = mt_name_conversion(pt(p).atropos.names,name);
            datlas = pt(p).dkt.label;
            dlabels = mt_name_conversion(pt(p).dkt.names,name);
            assert(isequal(alabels,dlabels))
    
            % convert musc to hup
            alabels = convert_musc_labels_to_hup(alabels);
            dlabels = convert_musc_labels_to_hup(dlabels);
        
            % Fix for HUP162
            if strcmp(pt(p).name,'HUP162')
                alabels = [alabels;'LA12';'LB12';'LC12';'RA12';'RB12';'RC12'];
                aatlas = [aatlas;'nan';'nan';'nan';'nan';'nan';'nan'];
                datlas = [datlas;'nan';'nan';'nan';'nan';'nan';'nan'];
            end
    
            % Fix for MP0037
            if strcmp(pt(p).name,'MP0037')
                alabels = [alabels;'RC5'];
                aatlas = [aatlas;'nan'];
                datlas = [datlas;'nan'];
            end
    
    
            [Lia,~] = ismember(alabels,clean_labels);
            [~,Locb] = ismember(clean_labels,alabels(Lia));
            txf_fcn1 = @(x) x(Lia);
            txf_fcn2 = @(x) x(Locb);
            txf_fcn = @(x) txf_fcn2(txf_fcn1(x));
    
            assert(isequal(clean_labels,txf_fcn(alabels)))
            aatlas = txf_fcn(aatlas);
            datlas = txf_fcn(datlas);
            
            if 0
                table(clean_labels,aatlas,datlas)
            end
            %}
    
            
        end
    end

    % find disconnected periods
    bp_machine_bb = squeeze(bp(:,1,:,6)); % get machine reference broadband bandpower
    bp_machine_bb = squeeze(nansum(abs(bp_machine_bb),2))/size(bp,3);% sum abs across channels, divide by nchs
    low_amplitude = bp_machine_bb < 1e-10; % these periods are likely disconnected
    all_disconnected(p,:) = low_amplitude;

    % also load the sleep stages
    stage = load([sleep_stage_path,name,'/sleep_stage.mat']);
    sout = stage.sout;

    % Get times of sleep transitions
    st_dates= sout.Summary(2:end,2);
    st_times = sout.Summary(2:end,3);
    stage = sout.Summary(2:end,4);

    % convert to seconds in ieeg
    seeg_secs = convert_transitions_to_ieeg_secs(st_dates,st_times);

    % get the sleep stages for appropriate times
    [~,seeg_stage] = match_seeg_times(times(:,1),ones(size(times,1),1),seeg_secs,stage);

    sleep = strcmp(seeg_stage,'N3') | strcmp(seeg_stage,'N2');% | ...
        %strcmp(seeg_stage,'N1') | strcmp(seeg_stage,'R');
    wake = strcmp(seeg_stage,'W');

    all_n_wake_sleep(p,:) = [sum(wake) sum(sleep)];
    all_n_wake_sleep_connected(p,:) = [sum(wake&~low_amplitude) sum(sleep&~low_amplitude)];
    

    % Save all, wake, and sleep
    for im = 1:3
        all_coh{p,im,1} = squeeze(nanmean(coh(~low_amplitude,im,:,:,:),1));
        all_coh{p,im,2} = squeeze(nanmean(coh(wake&~low_amplitude,im,:,:,:),1));
        all_coh{p,im,3} = squeeze(nanmean(coh(sleep&~low_amplitude,im,:,:,:),1));

        all_plv{p,im,1} = squeeze(nanmean(plv(~low_amplitude,im,:,:,:),1));
        all_plv{p,im,2} = squeeze(nanmean(plv(wake&~low_amplitude,im,:,:,:),1));
        all_plv{p,im,3} = squeeze(nanmean(plv(sleep&~low_amplitude,im,:,:,:),1));

        all_re{p,im,1} = squeeze(nanmean(re(~low_amplitude,im,:,:,:),1));
        all_re{p,im,2} = squeeze(nanmean(re(wake&~low_amplitude,im,:,:,:),1));
        all_re{p,im,3} = squeeze(nanmean(re(sleep&~low_amplitude,im,:,:,:),1));
    
        all_spikes{p,im,1} = squeeze(nanmean(spike_counts(~low_amplitude,im,:),1));
        all_spikes{p,im,2} = squeeze(nanmean(spike_counts(wake&~low_amplitude,im,:),1));
        all_spikes{p,im,3} = squeeze(nanmean(spike_counts(sleep&~low_amplitude,im,:),1));

    
        all_pearson{p,im,1} = squeeze(nanmean(pc(~low_amplitude,im,:,:),1));
        all_pearson{p,im,2} = squeeze(nanmean(pc(wake&~low_amplitude,im,:,:),1));
        all_pearson{p,im,3} = squeeze(nanmean(pc(sleep&~low_amplitude,im,:,:),1));

        all_pearson_squared{p,im,1} = squeeze(nanmean(pc_squared(~low_amplitude,im,:,:),1));
        all_pearson_squared{p,im,2} = squeeze(nanmean(pc_squared(wake&~low_amplitude,im,:,:),1));
        all_pearson_squared{p,im,3} = squeeze(nanmean(pc_squared(sleep&~low_amplitude,im,:,:),1));

        all_xcor{p,im,1} = squeeze(nanmean(xcor(~low_amplitude,im,:,:),1));
        all_xcor{p,im,2} = squeeze(nanmean(xcor(wake&~low_amplitude,im,:,:),1));
        all_xcor{p,im,3} = squeeze(nanmean(xcor(sleep&~low_amplitude,im,:,:),1));
    
        all_lags{p,im,1} = squeeze(nanmean(lags(~low_amplitude,im,:,:),1));
        all_lags{p,im,2} = squeeze(nanmean(lags(wake&~low_amplitude,im,:,:),1));
        all_lags{p,im,3} = squeeze(nanmean(lags(sleep&~low_amplitude,im,:,:),1));

        all_rl{p,im,1} = squeeze(nanmean(rl(~low_amplitude,im,:),1));
        all_rl{p,im,2} = squeeze(nanmean(rl(wake&~low_amplitude,im,:),1));
        all_rl{p,im,3} = squeeze(nanmean(rl(sleep&~low_amplitude,im,:),1));

        all_se{p,im,1} = squeeze(nanmean(se(~low_amplitude,im,:),1));
        all_se{p,im,2} = squeeze(nanmean(se(wake&~low_amplitude,im,:),1));
        all_se{p,im,3} = squeeze(nanmean(se(sleep&~low_amplitude,im,:),1));
    
        all_bp{p,im,1} = squeeze(nanmean(bp(~low_amplitude,im,:,:),1));
        all_bp{p,im,2} = squeeze(nanmean(bp(wake&~low_amplitude,im,:,:),1));
        all_bp{p,im,3} = squeeze(nanmean(bp(sleep&~low_amplitude,im,:,:),1));

        all_rel_bp{p,im,1} = squeeze(nanmean(rel_bp(~low_amplitude,im,:,:),1));
        all_rel_bp{p,im,2} = squeeze(nanmean(rel_bp(wake&~low_amplitude,im,:,:),1));
        all_rel_bp{p,im,3} = squeeze(nanmean(rel_bp(sleep&~low_amplitude,im,:,:),1));

        all_ll{p,im,1} = squeeze(nanmean(ll(~low_amplitude,im,:,:),1));
        all_ll{p,im,2} = squeeze(nanmean(ll(wake&~low_amplitude,im,:,:),1));
        all_ll{p,im,3} = squeeze(nanmean(ll(sleep&~low_amplitude,im,:,:),1));

        %% Do spike subsampling analysis
        connected_all = (find(~low_amplitude))';
        connected_sleep = (find(~low_amplitude & sleep))';
        connected_wake = (find(~low_amplitude & wake))';

        % Loop over durations
        for id = 1:length(durations)

            % Method 1: random independent samples (with replacement)
            % If duration requested exceeds available data, it will just
            % resample the same thing many times. I think this is ok.
            rand_all = randi(length(connected_all),durations(id),samples_per_duration); 
            rand_all = connected_all(rand_all);

            if isempty(connected_sleep)
                rand_sleep = []; 
            else
                rand_sleep = randi(length(connected_sleep),durations(id),samples_per_duration); 
                rand_sleep = connected_sleep(rand_sleep);
            end
            

            if isempty(connected_wake)
                rand_wake = []; 
            else
                rand_wake = randi(length(connected_wake),durations(id),samples_per_duration); 
                rand_wake = connected_wake(rand_wake);
            end
            
            for is = 1:samples_per_duration
                
                
                spikes_subsample{p,im,1,1,id,is} = squeeze(nanmean(spike_counts(rand_all(:,is),im,:),1)); % average across the samples in the duration
                
                if ~isempty(rand_wake)
                    spikes_subsample{p,im,2,1,id,is} = squeeze(nanmean(spike_counts(rand_wake(:,is),im,:),1));
                end

                if ~isempty(rand_sleep)
                    spikes_subsample{p,im,3,1,id,is} = squeeze(nanmean(spike_counts(rand_sleep(:,is),im,:),1));
                end
            end

            % Method 2: continuous samples. If I request more than I have,
            % it will take what it can.
            if 0
                rand_all = randi(max(1,length(connected_all)-durations(id)+1),1,samples_per_duration);
                if isempty(connected_wake)
                    rand_wake = [];
                else
                    rand_wake = randi(max(length(connected_wake)-durations(id)+1,1),1,samples_per_duration);
                end
    
                if isempty(connected_sleep)
                    rand_sleep = [];
                else
                    rand_sleep = randi(max(length(connected_sleep)-durations(id)+1,1),1,samples_per_duration);
    
                end
    
                
                amount_to_add = repmat((1:durations(id))',1,samples_per_duration); % get a continuous length
    
                actual_times_all = rand_all+amount_to_add;
                actual_times_all = max(actual_times_all,1);
                actual_times_all = min(actual_times_all,length(connected_all));
                rand_all_added = connected_all(actual_times_all); % get the ordered connected periods
    
                if ~isempty(rand_wake)
                    actual_times_wake = rand_wake+amount_to_add;
                    actual_times_wake = max(actual_times_wake,1);
                    actual_times_wake = min(actual_times_wake,length(connected_wake));
                    rand_wake_added = connected_wake(actual_times_wake);
                else
                    rand_wake_added = [];
                end
    
                if ~isempty(rand_sleep)
    
                    actual_times_sleep = rand_sleep+amount_to_add;
                    actual_times_sleep = max(actual_times_sleep,1);
                    actual_times_sleep = min(actual_times_sleep,length(connected_sleep));
                    rand_sleep_added = connected_sleep(actual_times_sleep);
                else
                    rand_sleep_added = [];
                end
    
                
                
                
    
                for is = 1:samples_per_duration
                    spikes_subsample{p,im,1,2,id,is} = squeeze(nanmean(spike_counts(rand_all_added(:,is),im,:),1));
    
                    if ~isempty(rand_wake_added)
                        spikes_subsample{p,im,2,2,id,is} = squeeze(nanmean(spike_counts(rand_wake_added(:,is),im,:),1));
                    end
    
                    if ~isempty(rand_sleep_added)
                        spikes_subsample{p,im,3,2,id,is} = squeeze(nanmean(spike_counts(rand_sleep_added(:,is),im,:),1));
                    end
                end
            end

        end

        if 0
            % test the subsampling
            figure
            tiledlayout(1,2)
            for is = 1:2
                nexttile
                a = squeeze(spikes_subsample(p,im,1,is,:,:));
                for id = 1:length(durations)
                    b = cell2mat(a(id,:));
                    rho = corr(b,'rows','pairwise');
                    drho = logical(eye(size(rho)));
                    rho(drho) = nan;
                    median_rho = nanmedian(rho,'all');
                    plot(id,median_rho,'o','markersize',10)
                    hold on
                end
            end
        end

        %% Also measure STD

        all_coh_iqr{p,im,1} = squeeze(nanstd(coh(~low_amplitude,im,:,:,:),[],1));
        all_coh_iqr{p,im,2} = squeeze(nanstd(coh(wake&~low_amplitude,im,:,:,:),[],1));
        all_coh_iqr{p,im,3} = squeeze(nanstd(coh(sleep&~low_amplitude,im,:,:,:),[],1));

        all_plv_iqr{p,im,1} = squeeze(nanstd(plv(~low_amplitude,im,:,:,:),[],1));
        all_plv_iqr{p,im,2} = squeeze(nanstd(plv(wake&~low_amplitude,im,:,:,:),[],1));
        all_plv_iqr{p,im,3} = squeeze(nanstd(plv(sleep&~low_amplitude,im,:,:,:),[],1));

        all_re_iqr{p,im,1} = squeeze(nanstd(re(~low_amplitude,im,:,:,:),[],1));
        all_re_iqr{p,im,2} = squeeze(nanstd(re(wake&~low_amplitude,im,:,:,:),[],1));
        all_re_iqr{p,im,3} = squeeze(nanstd(re(sleep&~low_amplitude,im,:,:,:),[],1));
    
        all_spikes_iqr{p,im,1} = squeeze(nanstd(spike_counts(~low_amplitude,im,:),[],1));
        all_spikes_iqr{p,im,2} = squeeze(nanstd(spike_counts(wake&~low_amplitude,im,:),[],1));
        all_spikes_iqr{p,im,3} = squeeze(nanstd(spike_counts(sleep&~low_amplitude,im,:),[],1));
    
        all_pearson_iqr{p,im,1} = squeeze(nanstd(pc(~low_amplitude,im,:,:),[],1));
        all_pearson_iqr{p,im,2} = squeeze(nanstd(pc(wake&~low_amplitude,im,:,:),[],1));
        all_pearson_iqr{p,im,3} = squeeze(nanstd(pc(sleep&~low_amplitude,im,:,:),[],1));

        all_pearson_squared_iqr{p,im,1} = squeeze(nanstd(pc_squared(~low_amplitude,im,:,:),[],1));
        all_pearson_squared_iqr{p,im,2} = squeeze(nanstd(pc_squared(wake&~low_amplitude,im,:,:),[],1));
        all_pearson_squared_iqr{p,im,3} = squeeze(nanstd(pc_squared(sleep&~low_amplitude,im,:,:),[],1));

        all_xcor_iqr{p,im,1} = squeeze(nanstd(xcor(~low_amplitude,im,:,:),[],1));
        all_xcor_iqr{p,im,2} = squeeze(nanstd(xcor(wake&~low_amplitude,im,:,:),[],1));
        all_xcor_iqr{p,im,3} = squeeze(nanstd(xcor(sleep&~low_amplitude,im,:,:),[],1));
    
        all_lags_iqr{p,im,1} = squeeze(nanstd(lags(~low_amplitude,im,:,:),[],1));
        all_lags_iqr{p,im,2} = squeeze(nanstd(lags(wake&~low_amplitude,im,:,:),[],1));
        all_lags_iqr{p,im,3} = squeeze(nanstd(lags(sleep&~low_amplitude,im,:,:),[],1));

        all_rl_iqr{p,im,1} = squeeze(nanstd(rl(~low_amplitude,im,:),[],1));
        all_rl_iqr{p,im,2} = squeeze(nanstd(rl(wake&~low_amplitude,im,:),[],1));
        all_rl_iqr{p,im,3} = squeeze(nanstd(rl(sleep&~low_amplitude,im,:),[],1));

        all_se_iqr{p,im,1} = squeeze(nanstd(se(~low_amplitude,im,:),[],1));
        all_se_iqr{p,im,2} = squeeze(nanstd(se(wake&~low_amplitude,im,:),[],1));
        all_se_iqr{p,im,3} = squeeze(nanstd(se(sleep&~low_amplitude,im,:),[],1));
    
        all_bp_iqr{p,im,1} = squeeze(nanstd(bp(~low_amplitude,im,:,:),[],1));
        all_bp_iqr{p,im,2} = squeeze(nanstd(bp(wake&~low_amplitude,im,:,:),[],1));
        all_bp_iqr{p,im,3} = squeeze(nanstd(bp(sleep&~low_amplitude,im,:,:),[],1));

        all_rel_bp_iqr{p,im,1} = squeeze(nanstd(rel_bp(~low_amplitude,im,:,:),[],1));
        all_rel_bp_iqr{p,im,2} = squeeze(nanstd(rel_bp(wake&~low_amplitude,im,:,:),[],1));
        all_rel_bp_iqr{p,im,3} = squeeze(nanstd(rel_bp(sleep&~low_amplitude,im,:,:),[],1));

        all_ll_iqr{p,im,1} = squeeze(nanstd(ll(~low_amplitude,im,:,:),[],1));
        all_ll_iqr{p,im,2} = squeeze(nanstd(ll(wake&~low_amplitude,im,:,:),[],1));
        all_ll_iqr{p,im,3} = squeeze(nanstd(ll(sleep&~low_amplitude,im,:,:),[],1));
    end

    

end

%% Save
out.all_names = all_names;
out.all_labels = all_labels;
out.all_coh = all_coh;
out.all_spikes = all_spikes;
out.all_rl = all_rl;
out.all_pearson = all_pearson;
out.all_pearson_squared = all_pearson_squared;
out.all_bp = all_bp;
out.all_plv = all_plv;
out.all_re = all_re;
out.all_se = all_se;
out.all_rel_bp = all_rel_bp;
out.all_xcor = all_xcor;
out.all_lags = all_lags;
out.all_ll = all_ll;
out.all_coh_iqr = all_coh_iqr;
out.all_spikes_iqr = all_spikes_iqr;
out.all_rl_iqr = all_rl_iqr;
out.all_pearson_iqr = all_pearson_iqr;
out.all_pearson_squared_iqr = all_pearson_squared_iqr;
out.all_bp_iqr = all_bp_iqr;
out.all_plv_iqr = all_plv_iqr;
out.all_re_iqr = all_re_iqr;
out.all_se_iqr = all_se_iqr;
out.all_rel_bp_iqr = all_rel_bp_iqr;
out.all_xcor_iqr = all_xcor_iqr;
out.all_lags_iqr = all_lags_iqr;
out.all_ll_iqr = all_ll_iqr;
out.all_disconnected = all_disconnected;

out.all_surgery = all_surgery;
out.all_resec_lat = all_resec_lat;
out.all_resec_loc = all_resec_loc;
out.all_ablate_lat = all_ablate_lat;
out.all_ablate_loc = all_ablate_loc;
out.all_engel = all_engel;
out.all_ilae = all_ilae;
out.all_soz_loc = all_soz_loc;
out.all_soz_lat = all_soz_lat;
%{
out.all_atropos = all_atropos;
out.all_dkt = all_dkt;
%}
out.all_atropos = [];
out.all_dkt = [];

out.all_n_wake_sleep = all_n_wake_sleep;
out.all_n_wake_sleep_connected = all_n_wake_sleep_connected;
out.spikes_subsample = spikes_subsample;
out.durations = durations;

save([out_folder,'mt_out_epilepsy_laterality.mat'],'out')

end