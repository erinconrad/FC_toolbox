function plot_example_detections

%% General parameters
n_sp = 50;
surround = 1;
overwrite = 1;

%% Locations
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
spike_folder = [results_folder,'all_out/'];
data_folder = [locations.main_folder,'data/'];
out_folder = [results_folder,'sp_validation/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;


whichPts = [];
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


for p = whichPts
    pt_name = pt(p).name;
    
    
    for im = 1:2
        outname = [out_folder,sprintf('%s_montage%d',pt_name,im)];
        if exist(outname,'file') ~= 0
  
            if overwrite == 0
                fprintf('\nSkipping %s\n',pt_name);
                continue
            
            else
                fprintf('\nDoing %s\n',pt_name);
            end
        else
            fprintf('\nDoing %s\n',pt_name);
        end

        %% Load spike file
        out = load([spike_folder,sprintf('%s_pc.mat',pt_name)]);
        out = out.pc;

        %% concatenate all spikes into one long thing
        % Include an extra column for the file index and block
        all_spikes = [];
        for f = 1:length(out.file)

            for h = 1:length(out.file(f).run)
                gdf = out.file(f).run(h).data.montage(im).spikes;

                all_spikes = [all_spikes;gdf,...
                    repmat(f,size(gdf,1),1),...
                    repmat(h,size(gdf,1),1)];
            end
        end

        %% initialize figure
        figure
        set(gcf,'position',[0 0 1400 1000])
        tiledlayout(ceil(n_sp/5),5,'tilespacing','tight','padding','tight');

        % Loop over spikes
        for i = 1:n_sp

            %% Randomly pick spike
            sp = randi(size(all_spikes,1));

            %% Get info about the spike
            f = all_spikes(sp,3);
            h = all_spikes(sp,4);
            fs = out.file(f).run(h).data.fs;
            sp_index = all_spikes(sp,2);
            run_start = out.file(f).run(h).run_times(1);
            sp_time = (sp_index-1)/fs + run_start;        
            sp_ch = all_spikes(sp,1);


            fname = pt(p).ieeg.file(f).name;
            which_chs = find(out.file(f).run(h).data.montage(im).is_run);


            %% Get the EEG data
            run_times = [sp_time - surround,sp_time+surround];
            data = download_ieeg_data(fname, login_name, pwfile, run_times,1);
            values = data.values;
            chLabels = data.chLabels;
            sp_index = surround*fs;

            clean_labs = decompose_labels(chLabels);
            if im == 2
                [values,labels] = car_montage(values,which_chs,clean_labs);
            else
                [values,~,labels] =...
                bipolar_montage(values,chLabels,which_chs,[],[]);
            end

            %% Plot data
            nexttile
            plot(linspace(0,surround*2,size(values,1)),values(:,sp_ch),'linewidth',2);
            hold on
            plot(surround,values(round(sp_index),sp_ch),'o','markersize',10)
                title(sprintf('Spike %d %1.1f s %s file %d',...
                    sp,sp_time,labels{sp_ch},f),'fontsize',10)

            yticklabels([])
            set(gca,'fontsize',10)


        end

        print(outname,'-djpeg');
        close(gcf)
    end
    
end

end