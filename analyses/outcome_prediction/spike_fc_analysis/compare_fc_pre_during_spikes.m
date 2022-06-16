function compare_fc_pre_during_spikes(overwrite)

N = 100;
surround = 3;
tw = 2;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/data/spike_corrs/'];
spikes_folder = [results_folder,'all_out/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end
login_name = locations.ieeg_login;
pwfile = locations.ieeg_pw_file;
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));

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
good_pts = T.GoodCARSpikes;
good_pts = good_pts(~isnan(good_pts));
good_pt_names = T.Var14;
npts = length(good_pts);

for l = 1:npts
    j = good_pts(l);
    name = pt(j).name;
    rid = pt(j).rid;
    
    if exist([out_folder,name,'.mat'],'file') ~= 0 && overwrite == 0
        
        fprintf('\nSkipping %s\n',name);
        continue
    else
        fprintf('\nDoing %s\n',name);
    end
    
    
    %% Load the spike file
    fname = [spikes_folder,name,'_pc.mat'];
    if ~exist(fname,'file')
        fprintf('\nCannot find spike file for %s, skipping...\n',name);
        continue
    end
    %% Load spike file
    out = load(fname);
    out = out.pc;
    
    %% concatenate all spikes into one long thing
        % Include an extra column for the file index and block
    all_spikes = [];
    for f = 1:length(out.file)


        for h = 1:length(out.file(f).run)
            gdf = out.file(f).run(h).data.montage(2).spikes;

            all_spikes = [all_spikes;gdf,...
                repmat(h,size(gdf,1),1),...
                repmat(f,size(gdf,1),1)];
        end

    end
    
    %% Pick a random sample of these spikes
    nspikes = size(all_spikes,1);
    which_spikes = randsample(nspikes,N);
    
    all_chs_corr = nan(N,2);
    sp_chs_corr = nan(N,2);
    single_ch_corr = nan(N,2);
    
    %% Loop over these spikes ang grab data
    for is = 1:N
        if mod(is,10) == 1
            fprintf('\nDoing spike %d of %d\n',is,N);
        end
        s = which_spikes(is);
        f = all_spikes(s,4);
        h = all_spikes(s,3);
        sp_index = all_spikes(s,2);
        sp_ch = all_spikes(s,1);
        
        fs = out.file(f).run(h).data.fs;
        run_start = out.file(f).run(h).run_times(1);
        which_chs = find(out.file(f).run(h).data.montage(2).is_run);
        sp_time = (sp_index-1)/fs + run_start;
        fname = pt(j).ieeg.file(f).name;
        
        %% Get the EEG data
        run_times = [sp_time - surround,sp_time+surround];
        data = download_ieeg_data(fname, login_name, pwfile, run_times,1);
        values = data.values;
        chLabels = data.chLabels;
        sp_index = surround*fs;

        clean_labels = decompose_labels(chLabels,name);
        
        %% Find non-intracranial chs
        non_intracranial = find_non_intracranial(clean_labels);
        which_chs = find(~non_intracranial); % channels to do analysis on
        
        %% Reject bad channels
        [bad,details] = identify_bad_chs(values,which_chs,chLabels,fs);
        which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

        %% CAR
        [values,labels] = car_montage(values,which_chs,clean_labels);
        
         % make non run channels nans
        is_run = ismember((1:length(clean_labels))',which_chs);
        values(:,~is_run) = nan;

        % filters
        values = notch_filter(values,fs);
        values = bandpass_filter(values,fs);
        
        %% Spike detector
        gdf = detector_alt(values,fs);
        if isempty(gdf)
            continue
        end
        sp_chs = unique(gdf(:,1));
        
        %% Get networks
        pre_spike_period = 1:round(fs*tw); %3 seconds before spike to 1 second before spike
        spike_period = round(fs*tw) + 1 : round(fs*2*tw); % 1 second before spike to 1 second after spike
        
        
        
        % pre spike all channels
        clip = values(pre_spike_period,:);
        pc = corrcoef(clip);
        pc(logical(eye(size(pc)))) = nan;
        all_chs_corr(is,1) = nanmean(pc,'all');
        
        % spike all channels
        clip = values(spike_period,:);
        pc = corrcoef(clip);
        pc(logical(eye(size(pc)))) = nan;
        all_chs_corr(is,2) = nanmean(pc,'all');
        
        % pre spike spike channels
        clip = values(pre_spike_period,sp_chs);
        pc = corrcoef(clip);
        pc(logical(eye(size(pc)))) = nan;
        sp_chs_corr(is,1) = nanmean(pc,'all');
        
        % spike spike channels
        clip = values(spike_period,sp_chs);
        pc = corrcoef(clip);
        pc(logical(eye(size(pc)))) = nan;
        sp_chs_corr(is,2) = nanmean(pc,'all');
        
        % pre spike single spike channel
        clip = values(pre_spike_period,:);
        pc = corrcoef(clip);
        pc(logical(eye(size(pc)))) = nan;
        single_ch_corr(is,1) = nanmean(pc(:,sp_ch));
        
        % spike single spike channel
        clip = values(spike_period,:);
        pc = corrcoef(clip);
        pc(logical(eye(size(pc)))) = nan;
        single_ch_corr(is,2) = nanmean(pc(:,sp_ch));
        
        
        if 0
            figure
            nexttile
            plot(values(:,sp_ch))
            nexttile
            plot(values(spike_period,sp_ch))
            nexttile
            plot(values(pre_spike_period,sp_ch))
            
        end
        
    end
    
    spout.all_chs_corr = all_chs_corr;
    spout.sp_chs_corr = sp_chs_corr;
    spout.single_ch_corr = single_ch_corr;
    save([out_folder,name,'.mat'],'spout');

end