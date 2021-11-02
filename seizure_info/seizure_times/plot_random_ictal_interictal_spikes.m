function plot_random_ictal_interictal_spikes

%% Parameters
m = 2; % do not change
nspikes = 100;
spikes_per_fig = 20;
spikes_per_row = 2;
nfigs = nspikes/spikes_per_fig;
nrows = spikes_per_fig/spikes_per_row;
surround = 3;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'sz_removal/'];
spikes_folder = [results_folder,'all_out/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

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

% Create a cell array containing all spikes
spikes = cell(npts,2);
for i = 1:npts
    p = good_pts(i);
    name = pt(p).name;
    
    %% Load the spike file
    fname = [spikes_folder,name,'_pc.mat'];
    if ~exist(fname,'file')
        fprintf('\nCannot find spike file for %s, skipping...\n',name);
        continue
    end
    pc = load([spikes_folder,sprintf('%s_pc.mat',name)]);
    pc = pc.pc;
    
    interictal_spikes = [];
    ictal_spikes = [];
    
    nfiles = length(pc.file);
    
    for f = 1:nfiles
        
        if ~isfield(pt(p).ieeg.file(f),'sz_times')
            sz_times = [];
        else
            sz_times = pt(p).ieeg.file(f).sz_times;
        end
        
        nblocks = length(pc.file(f).run);
        
        for b = 1:nblocks
            gdf = (pc.file(f).run(b).data.montage(m).spikes);
            
            if isempty(gdf), continue; end
            if isempty(sz_times)
                interictal = logical(ones(size(gdf,1),1));
                ictal = ~interictal;
            else
                ictal = (any(gdf(:,2)' >= sz_times(:,1) & gdf(:,2)' <= sz_times(:,2),1))';
                interictal = ~ictal;
            end
            
            interictal_spikes = [interictal_spikes;...
                repmat([p,f,b],sum(interictal),1),gdf(interictal,:)];
            
            ictal_spikes = [ictal_spikes;...
                repmat([p,f,b],sum(ictal),1),gdf(ictal,:)];
                
            
        end
    end
    
    spikes{i,1} = interictal_spikes;
    spikes{i,2} = interictal_spikes;
    
end

%% Loop over nspikes
for i = 1:2 % ictal and interictal
    for i_f = 1:nfigs
        
        % initialize figure
        figure
        set(gcf,'position',[0 0 1400 1000])
        tiledlayout(nrows,spikes_per_row,'tilespacing','tight','padding','tight');
        
        % Loop over spikes in fig
        for s = 1:spikes_per_fig

            while 1
                % randomly pick a patient
                ip = randi(length(spikes));
                if size(spikes{ip,i},1) ~=0, break; end
            end

            % randomly pick a spike for that patient

            is = randi(size(spikes{ip,i},1));
            
            
            % get spike info
            p = spikes{ip,i}(is,1);
            f = spikes{ip,i}(is,2);
            b = spikes{ip,i}(is,3);
            sp_ch = spikes{ip,i}(is,4); 
            sp_index = spikes{ip,i}(is,5);
            name = pt(p).name;
            
            pc = load([spikes_folder,sprintf('%s_pc.mat',name)]);
            pc = pc.pc;
            fs = pt(p).ieeg.file(f).fs;
            which_chs = find(pc.file(f).run(b).data.montage(m).is_run);
            run_start = pc.file(f).run(b).run_times(1);
            sp_time = (sp_index-1)/fs + run_start;  
            
            fname = pt(p).ieeg.file(f).name;
            
            
            % Get the EEG data
            run_times = [sp_time - surround,sp_time+surround];
            data = download_ieeg_data(fname, login_name, pwfile, run_times,1);
            values = data.values;
            chLabels = data.chLabels;
            sp_index = surround*fs;
            
            clean_labs = decompose_labels(chLabels,name);
            if m == 2
                [values,labels] = car_montage(values,which_chs,clean_labs);
            else
                [values,~,labels] =...
                bipolar_montage(values,chLabels,which_chs,[],[],pt_name);
            end
            
            % filters
            values = notch_filter(values,fs);
            values = bandpass_filter(values,fs);
            
            % Plot data
            nexttile
            plot(linspace(0,surround*2,size(values,1)),values(:,sp_ch),'linewidth',2);
            hold on
            plot(surround,values(round(sp_index),sp_ch),'o','markersize',10)
                title(sprintf('%s %1.1f s %s file %d',...
                    name,sp_time,labels{sp_ch},f),'fontsize',10)

            yticklabels([])
            set(gca,'fontsize',10)
            
        end

        if i == 1
            ic_text = 'interictal';
        elseif i == 2
            ic_text = 'ictal';
        end
        
        outname = [ic_text,i_f];
        
        print([out_folder,outname],'-djpeg');
        close(gcf)
    end
end


end

