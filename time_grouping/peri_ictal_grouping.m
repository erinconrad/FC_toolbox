function out = peri_ictal_grouping(pc)

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% sz times file

nfiles = length(pc.file);

for f = 1:nfiles
    
    % loop over szs
    nszs = length(pc.file(f).sz);
    
    for s = 1:nszs
        
        nruns = length(pc.file(f).sz(s).run);
        
        % Prep matrices
        nmontages = length(pc.file(f).sz(s).run(1).data.montage);
        nchs = length(pc.file(f).sz(s).run(1).data(1).clean_labels);
        clean_labels = pc.file(f).sz(s).run(1).data(1).clean_labels;
        run_times = nan(nruns,2);
        
        for m = 1:nmontages
            net_montage{m} = nan(nchs*(nchs-1)/2,nruns);
            spikes_montage{m} = nan(nchs,nruns);
            ad_montage{m} = nan(nchs,nruns);
            coa_montage{m} = nan(nchs*(nchs-1)/2,nruns);
            rl_montage{m} = nan(nchs,nruns);
            %coi_montage{m} = nan(nchs,nruns);
            coi_global_montage{m} = nan(nruns,1);
        end
        
        for r = 1:nruns
        
            run_times(r,:) = (pc.file(f).sz(s).run(r).run_times);
            file_times(r) = mean(pc.file(f).sz(s).run(r).run_times);
            
            %% Get the data and calculate ns
            for m = 1:nmontages
                data = pc.file(f).sz(s).run(r).data.montage(m).net.data;
                fs = pc.file(f).sz(s).run(r).data.fs;
                is_run = pc.file(f).sz(s).run(r).data.montage(m).is_run;
                gdf = pc.file(f).sz(s).run(r).data.montage(m).spikes;
                ad = pc.file(f).sz(s).run(r).data.montage(m).ad;
                
                %% Network
                % unwrap
                data_uw = wrap_or_unwrap_adjacency_fc_toolbox(data);

                % make diagonal nans
                data_uw(logical(eye(size(data_uw)))) = nan;

                % Take R^2
                data_uw = (data_uw).^2;
                
                %% Spikes
                spikes = nan(nchs,1); % default to nan
                %coi_ch = nan(nchs,1);
                global_coi = nan;
                coa = nan(nchs,nchs);
                rl = nan(nchs,1);
                spikes(is_run) = 0; % default zero if we run it

                if ~isempty(gdf)
                    % Spike chs
                    spike_chs = gdf(:,1);

                    % Get counts
                    [C,~,ic]= unique(spike_chs);
                    a_counts = accumarray(ic,1);

                    spikes(C) = a_counts; % add counts

                    % Get spike coi
                    %[coi_ch,global_coi] = get_spike_coi(gdf,nchs,fs);
                    [coa,rl,global_coi] = build_sequences(gdf,nchs,fs);
                end
                
                % Fill up cell arrays
                net_montage{m}(:,r) = wrap_or_unwrap_adjacency_fc_toolbox(data_uw);
                coa_montage{m}(:,r) = wrap_or_unwrap_adjacency_fc_toolbox(coa);
                rl_montage{m}(:,r) = rl;
                spikes_montage{m}(:,r) = spikes;
                ad_montage{m}(:,r) = ad;
                %coi_montage{m}(:,r) = coi_ch;
                coi_global_montage{m}(r) = global_coi;
                
                
            end
            
        end
        
        % prep out
        for m = 1:nmontages
            out.file(f).sz(s).montage(m).net = net_montage{m};
            out.file(f).sz(s).montage(m).spikes = spikes_montage{m};
            out.file(f).sz(s).montage(m).ad = ad_montage{m};
            %out.file(f).montage(m).coi_ch = coi_montage{m};
            out.file(f).sz(s).montage(m).coi_global = coi_global_montage{m};
            out.file(f).sz(s).montage(m).labels = clean_labels;
            out.file(f).sz(s).montage(m).coa = coa_montage{m};
            out.file(f).sz(s).montage(m).rl = rl_montage{m};

        end
        out.file(f).sz(s).run_times = run_times;
        out.file(f).sz(s).clean_labels = clean_labels;
        
        
    end
    
end


end