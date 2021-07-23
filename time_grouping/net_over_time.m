function out = net_over_time(pc)


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
    nchs = length(pc.file(f).run(1).data(1).clean_labels);
    clean_labels = pc.file(f).run(1).data(1).clean_labels;
    run_center = nan(nruns,1);
    
    for m = 1:nmontages
        net_montage{m} = nan(nchs*(nchs-1)/2,nruns);
        ns_montage{m} = nan(nchs,nruns);
    end
    
    for r = 1:nruns
        
        run_center(r) = mean(pc.file(f).run(r).run_times);
        
        %% Get the data and calculate ns
        for m = 1:nmontages
            data = pc.file(f).run(r).data.montage(m).net.data;
            is_run = pc.file(f).run(r).data.montage(m).is_run;
            
            % unwrap
            data_uw = wrap_or_unwrap_adjacency_fc_toolbox(data);
            
            % make diagonal nans
            data_uw(logical(eye(size(data_uw)))) = nan;
            
            % Take R^2
            data_uw = (data_uw).^2;
            
            net_montage{m}(:,r) = wrap_or_unwrap_adjacency_fc_toolbox(data_uw);
            
            % get ns
            ns = nansum(data_uw,1);
            ns(sum(~isnan(data_uw),1)==0) = nan;
            ns_montage{m}(:,r) = ns;
            
        end
    end
    
    %% Remove dangling times (immediately adjacent to disconnected periods; look bad)
    for m = 1:nmontages
        data = net_montage{m};
        
        data = remove_dangling_times(data,nruns);
        
        net_montage{m} = data;
        
        ns = ns_montage{m};
        
        ns = remove_dangling_times(ns,nruns);
        ns_montage{m} = ns;
        
    end
    
    % prep out
    for m = 1:nmontages
        out.file(f).montage(m).net = net_montage{m};
        out.file(f).montage(m).ns = ns_montage{m};
        out.file(f).montage(m).labels = pc.file(f).run(1).data.montage(m).labels;
        
    end
    out.file(f).run_center = run_center;
    out.file(f).clean_labels = clean_labels;
    
    
    
end

end