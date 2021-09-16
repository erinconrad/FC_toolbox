function out = net_over_time(pc)

span_to_look = 3;
max_nans = 3;

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

nfiles = length(pc.file);
for f = 1:nfiles
    nruns = length(pc.file(f).run);
    file_name = pc.file(f).name;
    fs = pc.file(f).run(1).data.fs;
    
    % Prep matrices
    nmontages = length(pc.file(f).run(1).data.montage);
    nchs = length(pc.file(f).run(1).data(1).clean_labels);
    clean_labels = pc.file(f).run(1).data(1).clean_labels;
    run_center = nan(nruns,1);
    
    for m = 1:nmontages
        net_montage{m} = nan(nchs*(nchs-1)/2,nruns);
        spike_montage{m} = nan(nchs,nruns);
        ad_montage{m} = nan(nchs,nruns);
        coi_montage{m} = nan(nchs,nruns);
        coi_global_montage{m} = nan(nruns,1);
    end
    
    for r = 1:nruns
        
        run_center(r) = mean(pc.file(f).run(r).run_times);
        file_times(r) = mean(pc.file(f).run(r).run_times);
        
        %% Get the data and calculate ns
        for m = 1:nmontages
            data = pc.file(f).run(r).data.montage(m).net.data;
            is_run = pc.file(f).run(r).data.montage(m).is_run;
            if isfield(pc.file(f).run(r).data.montage(m),'spikes')
                gdf = pc.file(f).run(r).data.montage(m).spikes;
            else
                gdf = [];
            end
            if isfield(pc.file(f).run(r).data.montage(m),'ad')
                ad = pc.file(f).run(r).data.montage(m).ad;
            else
                ad = nan(nchs,1);
            end
            
            %% Network
            % unwrap
            data_uw = wrap_or_unwrap_adjacency_fc_toolbox(data);
            
            % make diagonal nans
            data_uw(logical(eye(size(data_uw)))) = nan;
            
            % Take R^2
            data_uw = (data_uw).^2;
            
            %% Spikes
            spikes = nan(nchs,1); % default to nan
            coi_ch = nan(nchs,1);
            global_coi = nan;
            spikes(is_run) = 0; % default zero if we run it
            
            if ~isempty(gdf)
                % Spike chs
                spike_chs = gdf(:,1);

                % Get counts
                [C,~,ic]= unique(spike_chs);
                a_counts = accumarray(ic,1);

                spikes(C) = a_counts; % add counts
                
                % Get spike coi
                [coi_ch,global_coi] = get_spike_coi(gdf,nchs,fs);
            end
            
            % Fill up cell arrays
            net_montage{m}(:,r) = wrap_or_unwrap_adjacency_fc_toolbox(data_uw);
            spike_montage{m}(:,r) = spikes;
            ad_montage{m}(:,r) = ad;
            coi_montage{m}(:,r) = coi_ch;
            coi_global_montage{m}(r) = global_coi;

            
        end
    end
    
    %% Remove dangling times (immediately adjacent to disconnected periods; look bad)
    for m = 1:nmontages
        data = net_montage{m};        
        [data,all_adj_bad] = remove_dangling_times(data,nruns);        
        net_montage{m} = data;
        
        spikes_montage{m}(:,all_adj_bad) = nan;
        ad_montage{m}(:,all_adj_bad) = nan;
        coi_montage{m}(:,all_adj_bad) = nan;
        coi_global_montage{m}(all_adj_bad) = nan;
        %{
        spikes = spike_montage{m};
        spikes = remove_dangling_times(spikes,nruns);
        spike_montage{m} = spikes;
       
        ad = ad_montage{m};
        ad = remove_dangling_times(ad,nruns);
        ad_montage{m} = ad;
        %}
         
    end
    
    %% Also remove any surrounded by nans
    for m = 1:nmontages
        data = net_montage{m};
        data = remove_if_nans_surround(data,span_to_look,max_nans);
        net_montage{m} = data;
        
        
        spikes = spike_montage{m};
        spikes = remove_if_nans_surround(spikes,span_to_look,max_nans);
        spikes_montage{m} = spikes;
        
        ad = ad_montage{m};
        ad = remove_if_nans_surround(ad,span_to_look,max_nans);
        ad_montage{m} = ad;
        %}
        
    end
    
    % prep out
    for m = 1:nmontages
        out.file(f).montage(m).net = net_montage{m};
        out.file(f).montage(m).spikes = spikes_montage{m};
        out.file(f).montage(m).ad = ad_montage{m};
        out.file(f).montage(m).coi_ch = coi_montage{m};
        out.file(f).montage(m).coi_global = coi_global_montage{m};
        out.file(f).montage(m).labels = pc.file(f).run(1).data.montage(m).labels;
        
    end
    out.file(f).run_center = run_center;
    out.file(f).clean_labels = clean_labels;
    
    
    
end

end