function out = individual_run_mt(file_path,pt,which_pt,meta,which_edf_file,overlap_log_file,mT,attempt_remove_oob)

%% Parameters
show_data = 0;
test_flag = 0;
tw = 2; % time window for network calculations.


%% Load the edf file
C = strsplit(file_path,'/');
name = C{end-1};
info = edfinfo(file_path);

%% Get basic info
samples_per_time = info.NumSamples(1,1);
num_samples = samples_per_time*info.NumDataRecords;
fs = round(info.NumSamples(1,1)/seconds(info.DataRecordDuration));
labels = cellstr(info.SignalLabels);

name = pt(which_pt).name;

%% get allowable electrodes (only LA, LB, LC, RA, RB, RC electrodes)
potentially_allowable_labels = get_allowable_elecs(name);

%% Find which of these electrodes are outside the brain
if isempty(pt(which_pt).dkt)

else
    % Get atlas labels
    dkt = pt(which_pt).dkt.label;
    atropos = pt(which_pt).atropos.label;
    atlas_elec_names = pt(which_pt).atropos.names;
    
    % convert atlas names to the A, B, C convention
    atlas_elec_names = mt_name_conversion(atlas_elec_names,name);
    
    if attempt_remove_oob
        outside_brain = zeros(length(potentially_allowable_labels),1);
        for i = 1:length(potentially_allowable_labels)
            % find the matching atlas elec name
            match = strcmp(potentially_allowable_labels{i},atlas_elec_names);
        
            if match == 0, continue; end
        
            if strcmp(dkt{match},'EmptyLabel') && (strcmp(atropos{match},'CSF') ...
                    || strcmp(atropos{match},'EmptyLabel'))
                outside_brain(i) = 1;
            end
        end
        outside_brain = logical(outside_brain);
        
        if 0
            table(pt(which_pt).atropos.names,atlas_elec_names,atropos,dkt)
        end
    
        % remove those outside brain
        potentially_allowable_labels(outside_brain) = [];
    end
end

%% Remove some potentially allowable labels if they aren't really targeting mesial temporal region for that patient
exr = strcmp(mT.name,name); assert(sum(exr==1));
exc = mT.exclude{exr};
if ~isempty(exc)
    C = strsplit(exc,', ');
    rm_allow = zeros(length(potentially_allowable_labels),1);
    for i = 1:length(C) % loop over exclusion labels
        rm_allow(contains(potentially_allowable_labels,C{i},'ignorecase',true)) = 1;
    end
    potentially_allowable_labels(rm_allow==1) = [];
end


%% Turn non-A, B, C MT electrodes into A, B, C
old_labels = labels;
labels = mt_name_conversion(labels,name);


%% Find labels that match allowable electrodes and have symmetric coverage
[allowed_labels,final_allowed_idx] = find_mt_symmetric_coverage(labels,potentially_allowable_labels);
if isempty(allowed_labels)
    out = [];
    return
end
nallowed = length(allowed_labels);
old_allowed = old_labels(final_allowed_idx); % restrict the old labels to be the acceptable ones in the same order

% confirm that the allowed labels and old allowed labels have matching
% order
assert(isequal(mt_name_conversion(old_allowed,name),allowed_labels))

%% Get the seizure times from the corresponding ieeg file
% Find which is the correct ieeg file
ieeg_file_name = meta.file_name;
found_it = 0;
for i = 1:length(pt(which_pt).ieeg.file)
    if strcmp(ieeg_file_name,pt(which_pt).ieeg.file(i).name)
        found_it = 1;
        ieeg_file = i;
        break
    end
end
assert(found_it==1)

% Get seizure times
if ~isfield(pt(which_pt).ieeg.file(ieeg_file),'sz_times')
    sz_times = [];
else
    sz_times = pt(which_pt).ieeg.file(ieeg_file).sz_times;
end

%% Re-align seizure times so that the start of this edf file is 0 
% Get edf file start time
edf_file_start_time = meta.times(which_edf_file,1);
sz_times_aligned = sz_times - edf_file_start_time;


%% Error if not left and right
if sum(contains(allowed_labels,'L'))==0 || sum(contains(allowed_labels,'R'))==0
    error('what')
end
%}

%% Initialize values
values = nan(num_samples,nallowed);

% Separately call edfread for each signal
for is = 1:nallowed
    curr_signal = old_allowed{is}; % pull the old label
    
    % Get timetable for that signal
    T = edfread(file_path,'SelectedSignals',curr_signal);
      
    % Get variable name (may be changed once put into table)
    Var1 = T.Properties.VariableNames{1};
    
    %% Convert the time table to an array
    temp_values = nan(num_samples,1);

    % Loop over segments
    for s = 1:size(T,1)
        seg = T.(Var1){s};

        % Where are we in the temp values
        start_idx = (s-1)*samples_per_time+1;
        end_idx = s*samples_per_time;

        % Fill up values
        temp_values(start_idx:end_idx) = seg;
    end
    
    %% Fill up values
    values(:,is) = temp_values;
    
end

%% Downsample to 256 hz
if abs(fs-256)>1
    %old_values = values; old_fs = fs; old_times = linspace(0,size(old_values,1)/old_fs,size(old_values,1));
    p = 256; q = fs;
    values = resample(values,p,q);
    num_samples = size(values,1);
    fs = 256;
end

%% Get times
times = linspace(0,num_samples/fs,num_samples);

%% Find the times
max_start = length(times) - fs*60-1; % must start 60 seconds before the end

%{
I tested the sz code looking at pt HUP170 edf file 22 sz time 120539.93. I
confirmed that the data for LA1 looks approximately the same on ieeg.org as
in the edf file. I also tried doing some random runs for this patient and
saw it did indeed sometimes produce overlap and so try again

if 0
    figure
    test_times = [498.93 543.9];
    test_indices = round(test_times*fs);
    plot(values(test_indices(1):test_indices(2),1))
end

%}



% Wrap in a while loop
while_counts = 0;
while 1
    % Take a random one minute segment
    rand_start = randi(round(max_start));
    rand_end = rand_start + round(fs*60);
    
    % Narrow values down to this time
    curr_times = times(rand_start:rand_end);
    curr_values = values(rand_start:rand_end,:);

    % See if these times overlap with any seizure times.
    start_1 = curr_times(1); end_1 = curr_times(end);
    overlap = 0;
    overlapping_sz = 0;
    overlapping_sz_end = 0;
    overlapping_sz_aligned = 0;
    overlapping_sz_aligned_end = 0;
    for s = 1:size(sz_times_aligned,1)
        if do_times_overlap(start_1,end_1,...
                sz_times_aligned(s,1),sz_times_aligned(s,2)) == 1
            overlap = 1;
            overlapping_sz = sz_times(s,1);
            overlapping_sz_end = sz_times(s,2);
            overlapping_sz_aligned = sz_times_aligned(s,1);
            overlapping_sz_aligned_end = sz_times_aligned(s,2);
            break
        end
    end

    % If there is overlap, try again
    if overlap == 1
        while_counts = while_counts+1;

        % I should also print this to a log because it shouldn't happen
        % often
        T = readtable(overlap_log_file);
        T = [T;{pt(which_pt).name,meta.files{which_edf_file},overlapping_sz,...
            overlapping_sz_end,overlapping_sz_aligned,...
            overlapping_sz_aligned_end,start_1,end_1,while_counts}];
        writetable(T,overlap_log_file);

        if while_counts > 20 % give up and skip the run!
            out = [];
            return
        end
        continue;
    else
        break % if not, break out of loop
    end

end

%% Reject bad channels
which_chs = 1:nallowed;
[bad,details] = identify_bad_chs(curr_values,which_chs,allowed_labels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on
%}

%% Convert MUSC labels into HUP labels for simplicity
allowed_labels = convert_musc_labels_to_hup(allowed_labels);

%% CAR reference
[car_values,car_labels] = car_montage(curr_values,which_chs,allowed_labels);
is_run_car = ismember((1:length(car_labels))',which_chs);

%% Machine reference
machine_values = curr_values;
machine_labels = cellfun(@(x) sprintf('%s-Ref',x),allowed_labels,'uniformoutput',false);
is_run_machine = ismember((1:length(car_labels))',which_chs);

%% Bipolar reference
[bipolar_values,~,bipolar_labels,chs_in_bipolar] = ...
    bipolar_montage_fc(curr_values,allowed_labels,[],[],name);
bad_bipolar = any(ismember(chs_in_bipolar,bad),2);
empty = cellfun(@(x) strcmp(x,'-'),bipolar_labels);
which_chs_bipolar = 1:size(chs_in_bipolar,1);
which_chs_bipolar(bad_bipolar|empty) = [];
is_run_bipolar = ismember((1:length(allowed_labels))',which_chs_bipolar);

%% Get inter-electrode distance matrix
dm = pseudo_distance_mt(allowed_labels);

%% prep freqs
freqs = get_frequencies; % confirmed last is broadband

%% Calculate network
% Loop over montages
for im = 1:3 
   
    if im == 1
        montage = 'machine';
        values = machine_values;
        curr_labels = machine_labels;
        is_run = is_run_machine;
    elseif im == 2
        montage = 'car';
        values = car_values;
        is_run = is_run_car;
        curr_labels = car_labels;
    elseif im == 3
        montage = 'bipolar';
        values = bipolar_values;
        is_run = is_run_bipolar;
        curr_labels = bipolar_labels;
    end

    % notch filter
    values = notch_filter(values,fs);

    % bandpass filter 0.5-80
    broadband = freqs(end,:);
    values = bandpass_any(values,fs,broadband,4);
    
    % make non run channels nans
    run_values = values;
    run_values(:,~is_run) = nan;
    skip = find(~is_run);

    % Turn nans within channel into mean of channel
    for ich = 1:size(run_values,2)
        run_values(isnan(run_values(:,ich)),ich) = nanmean(run_values(:,ich));
    end

    % line length
    ll = line_length(run_values);

    % cross correlation with 2 second time windows (very similar to full
    % thing)
    [xcor,lags] = cross_correlation(run_values,fs,tw,1);
    
    if test_flag
        % cross correlation with full window
        [xcor2,lags2] = cross_correlation(run_values,fs,tw,1); % Erin edited to add tw
        figure; set(gcf,'position',[10 10 900 400])
        t = tiledlayout(1,2);
        nexttile; turn_nans_gray(xcor); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        nexttile; turn_nans_gray(lags); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        title(t,'Cross correlation (max and lags), full time windows')

        figure; set(gcf,'position',[10 10 900 400])
        t = tiledlayout(1,2);
        nexttile; turn_nans_gray(xcor2); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        nexttile; turn_nans_gray(lags2); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        title(t,'Cross correlation (max and lags), short time window')

        figure
        plot(xcor(:),xcor2(:),'o')
    end

    % Spectral entropy
    %{
    Should I do this across frequency bands???
    %}
    se = spectral_entropy(values,fs,tw,1); % Erin edited to add tw
    
    if test_flag
        se2 = spectral_entropy(values,fs,tw,1);
        figure
        nexttile
        plot(se,'o')
        xticks(1:length(se))
        xticklabels(curr_labels)
        title('Spectral entropy, full tw')
        ylim([0.5 0.8])

        nexttile
        plot(se2,'o')
        xticks(1:length(se))
        xticklabels(curr_labels)
        title('Spectral entropy, short tw')
        ylim([0.5 0.8])

        nexttile
        plot(se,se2,'o')
    end

    % Relative entropy - long and short time windows different
    re = relative_entropy(run_values,fs,tw,1); %broadband at the end
    if test_flag
        re2 = relative_entropy(run_values,fs,tw,0);
        figure; set(gcf,'position',[10 10 1400 400])
        t = tiledlayout(1,size(re,3));
        for i = 1:size(re,3)
            if sum(~(isnan(re(:,:,i))|isinf(re(:,:,i))),'all') == 0
                nexttile; continue;
            else
                nexttile; turn_nans_gray(re(:,:,i)); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
            end
        end
        title(t,'Relative entropy, short tw')

        figure; set(gcf,'position',[10 10 1400 400])
        t = tiledlayout(1,size(re,3));
        for i = 1:size(re,3)
            if sum(~(isnan(re2(:,:,i))|isinf(re2(:,:,i))),'all') == 0
                nexttile; continue;
            else
                nexttile; turn_nans_gray(re2(:,:,i)); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
            end
        end
        title(t,'Relative entropy, full tw')

        figure
        for i = 1:size(re,3)
            curr_re = re(:,:,i); curr_re2 = re2(:,:,i);
            nexttile; plot(curr_re(:),curr_re2(:),'o')
        end
        
    end
    
    % PC
    [pc,pc_squared] = new_pearson_calc(run_values,fs,tw,1);
    if test_flag
        pc2 =  new_pearson_calc(run_values,fs,tw,0);
        figure
        turn_nans_gray(pc); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        title('Pearson correlation, full window')

        figure
        turn_nans_gray(pc2); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        title('Pearson correlation, short window')

        figure
        plot(pc(:),pc2(:),'o')
    end
    
    % Spikes
    gdf = detector_new_timing(run_values,fs); % Erin changed tmul from 19 to 17
    
    % Get alpha delta ratio
    ad_rat = calc_ad(run_values,fs);
    
    % Get bandpower (performs similarly for short and long tw.)
    [bp,rel_bp] = bp_calc_2(run_values,fs,[],tw,1); % broadband at the end
    if test_flag
        [bp2,rel_bp2] = bp_calc_2(run_values,fs,[],tw,1);
        figure; set(gcf,'position',[10 10 900 400])
        t = tiledlayout(1,2);
        nexttile
        turn_nans_gray(bp)
        yticks(1:size(xcor,1)); yticklabels(curr_labels);
        nexttile
        turn_nans_gray(rel_bp)
        yticks(1:size(xcor,1)); yticklabels(curr_labels);
        title(t,'Bandpower (absolute and relative), fulltw')

        figure; set(gcf,'position',[10 10 900 400])
        t = tiledlayout(1,2);
        nexttile
        turn_nans_gray(bp2)
        yticks(1:size(xcor,1)); yticklabels(curr_labels);
        nexttile
        turn_nans_gray(rel_bp2)
        yticks(1:size(xcor,1)); yticklabels(curr_labels);
        title(t,'Bandpower (absolute and relative), short tw')

        figure
        plot(bp(:,1),bp2(:,1),'o')
    end
    
    % Get coherence. Two second time window breaks it. 10 seconds about the same as 1 minute
    coh = faster_coherence_calc(run_values,fs,tw,1); % broadband at end
    if test_flag
        coh2 = faster_coherence_calc(run_values,fs,tw,0);
        figure; set(gcf,'position',[10 10 1400 400])
        t = tiledlayout(1,size(coh,3));
        for i = 1:size(coh,3)
        nexttile; turn_nans_gray(coh(:,:,i)); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        end
        title(t,'Coherence, full tw')

        figure; set(gcf,'position',[10 10 1400 400])
        t = tiledlayout(1,size(coh,3));
        for i = 1:size(coh,3)
        nexttile; turn_nans_gray(coh2(:,:,i)); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        end
        title(t,'Coherence, short tw')

        figure; set(gcf,'position',[10 10 1400 400])
        tiledlayout(1,size(coh,3));
        for i = 1:size(coh,3)
        nexttile; thing1 = coh(:,:,i); thing2 = coh2(:,:,i); plot(thing1,thing2,'o')
        end
    end

    % PLV - very similar short vs long time window.
    plv = plv_calc(run_values,fs,tw,1); % broadband at the end
    if test_flag
        plv2 = plv_calc(run_values,fs,tw,1);
        figure; set(gcf,'position',[10 10 1400 400])
        t = tiledlayout(1,size(plv,3));
        for i = 1:size(plv,3)
        nexttile; turn_nans_gray(plv(:,:,i)); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        end
        title(t,'PLV, full time window')

        figure; set(gcf,'position',[10 10 1400 400])
        t = tiledlayout(1,size(plv,3));
        for i = 1:size(plv,3)
        nexttile; turn_nans_gray(plv2(:,:,i)); yticks(1:size(xcor,1)); yticklabels(curr_labels); xticks(1:size(xcor,1)); xticklabels(curr_labels);
        end
        title(t,'PLV, short time window')

        figure
        tiledlayout(1,size(coh,3));
        for i = 1:size(plv,3)
        nexttile; thing1 = plv(:,:,i); thing2 = plv2(:,:,i); plot(thing1,thing2,'o')
        end

    end

    
    out.montage(im).name = montage;
    out.montage(im).bp = bp;
    out.montage(im).rel_bp = rel_bp;
    out.montage(im).pc = pc;
    out.montage(im).pc_squared = pc_squared;
    out.montage(im).plv = plv;
    out.montage(im).coh = coh;
    out.montage(im).gdf = gdf;
    out.montage(im).ad = ad_rat;
    out.montage(im).skip = skip;
    out.montage(im).is_run = is_run;
    out.montage(im).labels = curr_labels;
    out.montage(im).xcor = xcor;
    out.montage(im).lags = lags;
    out.montage(im).re = re;
    out.montage(im).ll = ll;
    out.montage(im).se = se;
    out.clean_labels = allowed_labels;
    out.fs = fs;
    out.times = [curr_times(1) curr_times(end)];
    out.idx = [rand_start rand_end];
    out.file_path = file_path;
    out.name = name;
    out.freqs = freqs;
    out.dm = dm;
    out.old_allowed_labels = old_allowed;

    if show_data
        tout.montage(im).values = values;
        tout.montage(im).name = montage;
    end

end


 if show_data
    ex_chs = [];
    only_run = 0;
    show_montage = 2;
    simple_plot(tout,out,ex_chs,show_montage,out.montage(show_montage).gdf,...
        only_run,out.montage(show_montage).skip)
    %pause
    %close(gcf)
    return

    clear tout

    m = 1;
    f = 1;
    curr_re = out.montage(m).re(:,:,f);
    curr_pc = out.montage(m).pc(:,:);
    curr_xcor = out.montage(m).xcor(:,:);
    curr_plv = out.montage(m).plv(:,:,f);

    figure; tiledlayout(1,3);
    for i = 1:3
        nexttile
        turn_nans_gray(out.montage(i).pc)
        title(out.montage(i).name)
        xticks(1:length(out.montage(i).labels))
        xticklabels(out.montage(i).labels)
        yticks(1:length(out.montage(i).labels))
        yticklabels(out.montage(i).labels)
    end

    figure
    tiledlayout(1,3)
    for i = 1:3
        nexttile
        turn_nans_gray(out.montage(i).plv(:,:,4))
        title(out.montage(i).name)
        xticks(1:length(out.montage(i).labels))
        xticklabels(out.montage(i).labels)
        yticks(1:length(out.montage(i).labels))
        yticklabels(out.montage(i).labels)
    end

    figure
    tiledlayout(1,3)
    for i = 1:3
        nexttile
        turn_nans_gray(out.montage(i).re(:,:,4))
        title(out.montage(i).name)
        xticks(1:length(out.montage(i).labels))
        xticklabels(out.montage(i).labels)
        yticks(1:length(out.montage(i).labels))
        yticklabels(out.montage(i).labels)
    end

    figure
    tiledlayout(1,3)
    for i = 1:3
        nexttile
        turn_nans_gray(out.montage(i).re_short_tw(:,:,4))
        title(out.montage(i).name)
        xticks(1:length(out.montage(i).labels))
        xticklabels(out.montage(i).labels)
        yticks(1:length(out.montage(i).labels))
        yticklabels(out.montage(i).labels)
    end
 end


end