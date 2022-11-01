function out = net_over_time(pc,pt,j)

%% Parameters
rm_ictal_spikes = 1; % remove spikes in seizures?
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
    
    % get coherence runs
    coherence_blocks = pt(j).ieeg.file(f).coherence_blocks;
    first_coh = find(coherence_blocks); first_coh = first_coh(1);
    
    % get nfreqs
    %if exist(
    nfreqs1 = 6;%size(pc.file(f).run(first_coh).cohere_out.montage(2).bp,2);
    nfreqs2 = 5;
    
    % Get seizure times
    if ~isfield(pt(j).ieeg.file(f),'sz_times')
        sz_times = [];
    else
        sz_times = pt(j).ieeg.file(f).sz_times;
        
    end

    if ~isfield(pt(j).ieeg.file(f),'sz_semiology')
        sz_semiology = {};
    else
        sz_semiology = pt(j).ieeg.file(f).sz_semiology;
    end
    
    % Get file start time
    file_start_time = pt(j).ieeg.file(f).start_time;
    
    % initialize cell arrays
    net_montage = cell(nmontages,1);
    spikes_montage = cell(nmontages,1);
    ad_montage = cell(nmontages,1);
    coa_montage = cell(nmontages,1);
    rl_montage = cell(nmontages,1);
    coi_global_montage = cell(nmontages,1);
    n_rm_ictal = zeros(nmontages,1);
    coh_montage = cell(nmontages,1);
    bp_montage = cell(nmontages,1);
    %seq = cell(nmontages,1);
    
    for m = 1:nmontages
        net_montage{m} = nan(nchs*(nchs-1)/2,nruns);
        coh_montage{m} = nan(nchs*(nchs-1)/2,nfreqs1,nruns);
        spikes_montage{m} = nan(nchs,nruns);
        ad_montage{m} = nan(nchs,nruns);
        bp_montage{m} = nan(nchs,nfreqs2,nruns);
        coa_montage{m} = nan(nchs*(nchs-1)/2,nruns);
        rl_montage{m} = nan(nchs,nruns);
        leader_montage{m} = nan(nchs,nruns);
        %coi_montage{m} = nan(nchs,nruns);
        coi_global_montage{m} = nan(nruns,1);
        seq_info{m} = nan(2,nruns);
        %seq{m} = [];
    end
    
    for r = 1:nruns
        
        run_center(r) = mean(pc.file(f).run(r).run_times);
        file_times(r) = mean(pc.file(f).run(r).run_times);
        run_start = pc.file(f).run(r).run_times(1);
        
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
            
            %% Get coherence and bandpower
            if coherence_blocks(r) == 1 &&  ...
                    isfield(pc.file(f).run(r),'cohere_out') && ~isempty(pc.file(f).run(r).cohere_out) % comment out last part when done
                coh = pc.file(f).run(r).cohere_out.montage(m).coh;
                bp = pc.file(f).run(r).cohere_out.montage(m).bp;
            else
                coh = nan(nchs,nchs,nfreqs1);
                bp = nan(nchs,nfreqs2);
            end
            
            
            
            %% Network
            % unwrap
            data_uw = wrap_or_unwrap_adjacency_fc_toolbox(data);
            
            % make diagonal nans
            data_uw(logical(eye(size(data_uw)))) = nan;
            
            % Take R^2
            %data_uw = (data_uw).^2;
            
            %% Default thigns to nan
            spikes = nan(nchs,1); % default to nan
            %coi_ch = nan(nchs,1);
            global_coi = nan;
            coa = nan(nchs,nchs);
            nseq_and_length = nan(2,1);
            rl = nan(nchs,1);
            leader = nan(nchs,1);
            
            %% Change defaults to zeros if we run it
            % this way I distinguish between zero spikes because we just
            % didn't run it from zero spikes because none were detected
            spikes(is_run) = 0; % default zero if we run it
            nseq_and_length = [0 0];
            coa = zeros(nchs,nchs);
            
            if rm_ictal_spikes && ~isempty(gdf)
                [gdf,n_ictal] = rm_spikes_in_seizures(gdf,fs,run_start,sz_times);
                n_rm_ictal(m) = n_rm_ictal(m) + n_ictal;
            end
            
            if ~isempty(gdf)
                % Spike chs
                spike_chs = gdf(:,1);

                % Get counts
                [C,~,ic]= unique(spike_chs);
                a_counts = accumarray(ic,1);

                spikes(C) = a_counts; % add counts
                
                % Get spike coi
                %[coi_ch,global_coi] = get_spike_coi(gdf,nchs,fs);
                [coa,rl,global_coi,seq_lengths,leader] = build_sequences(gdf,nchs,fs);
                nseq_and_length = [length(seq_lengths) mean(seq_lengths)];
            else
               
            end
            
            % Fill up cell arrays
            net_montage{m}(:,r) = wrap_or_unwrap_adjacency_fc_toolbox(data_uw);
            coa_montage{m}(:,r) = wrap_or_unwrap_adjacency_fc_toolbox(coa);
            coh_montage{m}(:,:,r) = wrap_or_unwrap_adjacency_fc_toolbox(coh);
            bp_montage{m}(:,:,r) = bp;
            %seq{m} = [seq{m};seq_matrix];
            rl_montage{m}(:,r) = rl;
            spikes_montage{m}(:,r) = spikes;
            ad_montage{m}(:,r) = ad;
            %coi_montage{m}(:,r) = coi_ch;
            coi_global_montage{m}(r) = global_coi;
            seq_info{m}(:,r) = nseq_and_length;
            leader_montage{m}(:,r) = leader;

            
        end
    end
    
    %% Remove dangling times (immediately adjacent to disconnected periods; look bad)
    for m = 1:nmontages
        data = net_montage{m};        
        [data,all_adj_bad] = remove_dangling_times(data,nruns);        
        net_montage{m} = data;
        
        spikes_montage{m}(:,all_adj_bad) = nan;
        ad_montage{m}(:,all_adj_bad) = nan;
        rl_montage{m}(all_adj_bad) = nan;
        coa_montage{m}(:,all_adj_bad) = nan;
        %coi_montage{m}(:,all_adj_bad) = nan;
        coi_global_montage{m}(all_adj_bad) = nan;
        seq_info{m}(:,all_adj_bad) = nan;
        leader_montage{m}(:,all_adj_bad) = nan;
        coh_montage{m}(:,:,all_adj_bad) = nan;
        bp_montage{m}(:,:,all_adj_bad) = nan;
        
         
    end
    
    %% Also remove any surrounded by nans
    for m = 1:nmontages
        data = net_montage{m};
        data = remove_if_nans_surround(data,span_to_look,max_nans);
        net_montage{m} = data;
        
        
        spikes = spikes_montage{m};
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
        %out.file(f).montage(m).coi_ch = coi_montage{m};
        out.file(f).montage(m).coi_global = coi_global_montage{m};
        out.file(f).montage(m).labels = pc.file(f).run(1).data.montage(m).labels;
        out.file(f).montage(m).coa = coa_montage{m};
        out.file(f).montage(m).rl = rl_montage{m};
        out.file(f).montage(m).n_rm_ictal = n_rm_ictal(m);
        out.file(f).montage(m).seq_info = seq_info{m};
        %out.file(f).montage(m).seq = seq{m};
        out.file(f).montage(m).leader_montage = leader_montage{m};
        out.file(f).montage(m).coh = coh_montage{m};
        out.file(f).montage(m).bp = bp_montage{m};
        
    end
    out.file(f).run_center = run_center;
    out.file(f).clean_labels = clean_labels;
    out.file(f).sz_times = sz_times;
    out.file(f).sz_semiology = sz_semiology;
    out.file(f).file_start_time = file_start_time;
    
    
    
end

end