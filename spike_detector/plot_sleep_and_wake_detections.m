function plot_sleep_and_wake_detections

%% General parameters
n_sp = 50;
surround = 1;

%% Locations
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
data_folder = [locations.main_folder,'data/'];
int_folder = [results_folder,'analysis/intermediate/'];
spike_folder = [results_folder,'all_out/'];
out_folder = [results_folder,'sleep_wake_validation/'];
scripts_folder = locations.script_folder;


validation_file = [scripts_folder,'spike_detector/Manual validation.xlsx'];
%% Get the indices of the patients with good spikes
T = readtable(validation_file);
good_pts = T.Var13;
good_pts = good_pts(~isnan(good_pts));
npts = length(good_pts);

if ~exist(out_folder,'dir'), mkdir(out_folder); end

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%{
if isempty(whichPts)
listing = dir([spike_folder,'*.mat']);
for i = 1:length(listing)
    C = listing(i).name;
    temp_name = strsplit(C,'_');
    temp_name = temp_name{1};
    for j = 1:length(pt)
        pt_name = pt(j).name;
        if strcmp(temp_name,pt_name)
            whichPts = [whichPts,j];
            break
        end
    end
end
end
%}

%% Do AD validation in order to get ADR to discriminate wake from sleep
roc_out = ad_validation;
disc = roc_out.disc;

for l = 1:npts
    p = good_pts(l);
    pt_name = pt(p).name;
    
    
    for im = 2
        
        
        %% Load spike file
        out = load([spike_folder,sprintf('%s_pc.mat',pt_name)]);
        out = out.pc;
        
        %% Load the corresponding summ file
        summ = load([int_folder,pt_name]);
        summ = summ.summ;
        
        %% Get the normalized ADR and file 
        ad = summ.ad;
        labels = summ.labels;
        ekg = find_non_intracranial(labels);
        ad = ad(~ekg,:);
        ad = nanmean(ad,1);
        fidx = summ.file_index;
        [sleep,wake] = find_sleep_wake(ad,[],disc);
        
        %% Figure out block for each time
        bidx = ones(length(fidx),1);
        for it = 2:length(bidx)
            
            % if it's the same file, increase the index
            if fidx(it) == fidx(it-1)
                bidx(it) = bidx(it-1)+1;
            else
                % if different, reset to 1
                bidx(it) = 1;
            end
        end
        
        
        
        %% Skip incomplete pts
        
        % Get corresponding pt
        for j = 1:length(pt)
            if strcmp(pt(j).name,pt_name)
                break
            end
        end
        
        %{
        % Skip the patient if it's incomplete
        if length(out.file) < length(pt(j).ieeg.file) || ...
                length(out.file(end).run) < size(pt(j).ieeg.file(end).run_times,1)
            fprintf('\n%s incomplete, skipping\n',pt_name);
            continue
        end
        %}
        

        %% concatenate all spikes into one long thing
        % Include an extra column for the file index and block
        all_spikes = [];
        for f = 1:length(out.file)

            
            for h = 1:length(out.file(f).run)
                gdf = out.file(f).run(h).data.montage(im).spikes;

                all_spikes = [all_spikes;gdf,...
                    repmat(h,size(gdf,1),1),...
                    repmat(f,size(gdf,1),1)];
            end
            
        end
        
        %% Add sleep/wake status to spike info
        nspikes = size(all_spikes,1);
        spikes_sw = nan(nspikes,3);
        for i = 1:nspikes
            
            % get current time index
            f = all_spikes(i,4);
            b = all_spikes(i,3);
            idx = conv_fb_to_idx(f,b,fidx,bidx);
            spikes_sw(i,:) = [sleep(idx), wake(idx), idx];
        end
        
        
        % add index
        sp_idx = (1:nspikes)';
        all_spikes = [all_spikes,sp_idx];
        
        sleep_spikes = all_spikes(spikes_sw(:,1)==1,:);
        wake_spikes = all_spikes(spikes_sw(:,2)==1,:);
        
        %% Some checks
        % confirm that sleep and wake sets are disjoint
        assert(isempty(intersect(sleep_spikes(:,end),wake_spikes(:,end))))
        
        % plot fidx and bidx
        if 0
            figure
            tiledlayout(2,1)
            nexttile
            plot(fidx)
            nexttile
            plot(bidx)
        end
        
       % plot ADRs in sleep and wake
       %{
       I tested this for HUP 100, which has two files, and it pretty
       convincingly separated the wake labeled spikes from the sleep
       labeled spikes
       %}
       if 0
           figure
           nsleep = size(sleep_spikes,1);
           sleep_idx = spikes_sw(spikes_sw(:,1)==1,3);
           plot(randn(nsleep,1)*0.05+1,ad(sleep_idx),'o');
           hold on
           nwake = size(wake_spikes,1);
           wake_idx = spikes_sw(spikes_sw(:,2)==1,3);
           plot(randn(nwake,1)*0.05+2,ad(wake_idx),'o');
           
       end

        %% initialize figure
        for is = 1:2
            figure
            set(gcf,'position',[0 0 1400 1000])
            tiledlayout(ceil(n_sp/5),5,'tilespacing','tight','padding','tight');

            if is == 1
                curr_spikes = sleep_spikes;
                sleep_text = 'Sleep';
            else
                curr_spikes = wake_spikes;
                sleep_text = 'Wake';
            end
            
            % Loop over spikes
            for i = 1:n_sp

                %% Randomly pick spike
                
                sp = randi(size(curr_spikes,1));

                %% Get info about the spike
                f = curr_spikes(sp,4);
                h = curr_spikes(sp,3);

                fs = out.file(f).run(h).data.fs;
                run_start = out.file(f).run(h).run_times(1);
                which_chs = find(out.file(f).run(h).data.montage(im).is_run);
                
                sp_index = curr_spikes(sp,2);

                sp_time = (sp_index-1)/fs + run_start;        
                sp_ch = curr_spikes(sp,1);
                fname = pt(p).ieeg.file(f).name;



                %% Get the EEG data
                run_times = [sp_time - surround,sp_time+surround];
                data = download_ieeg_data(fname, login_name, pwfile, run_times,1);
                values = data.values;
                chLabels = data.chLabels;
                sp_index = surround*fs;

                clean_labs = decompose_labels(chLabels,pt_name);
                if im == 2
                    [values,labels] = car_montage(values,which_chs,clean_labs);
                else
                    [values,~,labels] =...
                    bipolar_montage(values,chLabels,which_chs,[],[],pt_name);
                end

                % filters
                values = notch_filter(values,fs);
                values = bandpass_filter(values,fs);

                %% Plot data
                nexttile
                plot(linspace(0,surround*2,size(values,1)),values(:,sp_ch),'linewidth',2);
                hold on
                plot(surround,values(round(sp_index),sp_ch),'o','markersize',10)
                    title(sprintf('%s spike %d %1.1f s %s file %d',...
                        sleep_text,sp,sp_time,labels{sp_ch},f),'fontsize',10)

                yticklabels([])
                set(gca,'fontsize',10)


            end
            outname = [out_folder,sprintf('%s_montage%d_%s.jpg',pt_name,im,sleep_text)];
            print(outname,'-djpeg');
            close(gcf)
        end
    end
    
end

end


%% function to convert a file/block combo to single index
function idx = conv_fb_to_idx(f,b,fidx,bidx)

idx = 0;
for i = 1:f-1
    idx = idx + max(bidx(fidx==i)); % add the number of blocks from the preceding files
end

% also add current block idx
idx = idx + b;

end