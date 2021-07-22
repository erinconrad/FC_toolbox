function out = ns_over_time(pc)

%% Parameters
im = 1;
spacing = 20;
do_gui = 0;

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

nfiles = length(pc.file);
for f = 1:nfiles
    nruns = length(pc.file(f).run);
    file_name = pc.file(f).name;
    
    % Prep matrices
    nmontages = length(pc.file(f).run(1).data.montage);
    nchs = length(pc.file(1).run(1).data(1).clean_labels);
    clean_labels = pc.file(1).run(1).data(1).clean_labels;
    run_center = nan(nruns,1);
    
    for m = 1:nmontages
        
        data_montage{m} = nan(nchs,nruns);
        net_montage{m} = nan(nchs,nchs);
        
    end
    
    for r = 1:nruns
        
        run_center(r) = mean(pc.file(f).run(r).run_times);
        
        %% Get the data and calculate ns
        for m = 1:nmontages
            data = pc.file(f).run(r).data.montage(m).net.data;
            is_run = pc.file(f).run(r).data.montage(m).is_run;
            
            % Unwrap
            data = wrap_or_unwrap_adjacency_fc_toolbox(data);
            data = (data).^2;
            
            % Get ns
            ns = nansum(data,1);
            
            % Make ~is_run chs nans
            ns(~is_run) = nan;
            
            % Fill up cell array
            data_montage{m}(:,r) = ns;
            net_montage{m} = nanmean(cat(3,net_montage{m},data),3);
            
        end
    end
    
    %% Remove dangling times (immediately adjacent to disconnected periods; look bad)
    for m = 1:nmontages
        data = data_montage{m};
        
        data = remove_dangling_times(data,nruns);
        
        data_montage{m} = data;
        
    end
    
    % prep out
    for m = 1:nmontages
        out.file(f).montage(m).data = data_montage{m};
        out.file(f).montage(m).net = net_montage{m};
    end
    
    % raster plot
    if do_gui
        figure    
        set(gcf,'position',[10 10 1200 1000])
        turn_nans_gray(data_montage{im})
        xticks(1:spacing:nruns)
        xticklabels(run_center(1:spacing:nruns))
        yticks(1:nchs)
        yticklabels(clean_labels)

        while 1
            try
                [x,~] = ginput;
            catch
                break
            end
            curr_center = run_center(round(x));
            fprintf('\nShowing 15 seconds for %s start time %1.1fs\n',...
                file_name,curr_center-7.5);
            times = [curr_center - 7.5 curr_center + 7.5];
            which_net = 'pc';

            quick_run(file_name,times,which_net,0,3)

        end
    end
    
    
end

end