function [fc_atlas,spike_atlas,soz_atlas,coh_atlas] = ...
    rebuild_atlas(fc,spikes,atlas_labels,regions,atlas_nums,...
    labels,soz,coh,force_same_num,locs,lats)

%{
This function takes FC matrices, spikes, soz designations, along with a
    vector indicating what atlas region each electrode is in, and then it
    returns a new FC, spike, SOZ in atlas space
    
force_same_num is a 1 or 0 indicating whether to require the numbers of
electrodes contributing to each atlas region to be the same.
%}

npts = length(fc);
nregions = length(atlas_nums);
nfreqs = size(coh{1},3);

%% Remove ekg from atlas
for ip = 1:npts
    ekg = find_non_intracranial(atlas_labels{ip});
    atlas_labels{ip}(ekg) = [];
    regions{ip}(ekg) = [];
    
    assert(isequal(atlas_labels{ip},labels{ip})); % make sure electrode labels line up
end

% initialize stuff in atlas space
fc_atlas = nan(nregions,nregions,npts);
spike_atlas = nan(nregions,npts);
soz_atlas = zeros(nregions,npts);
coh_atlas = nan(nregions,nregions,nfreqs,npts);


% loop over patients
for ip = 1:npts
    
    curr_spikes = spikes{ip};
    curr_fc = fc{ip};
    curr_regions = regions{ip};
    curr_soz = soz{ip};
    curr_coh = coh{ip};
    
    %% Subsample channels
    subsample_regions = nan(length(regions{ip}),1);
    
    % Loop over regions
    for ir = 1:nregions
        % get this region num
        curr_num_i = atlas_nums(ir);
        
        % get the channels in this region
        curr_chs_i = ismember(curr_regions,curr_num_i);
        
        %% Get the channels in opposite region
        % current region name
        curr_loc_i = locs{ir};
        
        % opposite region
        opp_loc = find_opp_region(ir,curr_loc_i,locs,lats);
        
        % Get channels in opposite region
        curr_chs_i_opp = ismember(curr_regions,opp_loc);
        
        % Get smallest number of channels in i or opposite
        smallest_i = min([sum(curr_chs_i),sum(curr_chs_i_opp)]);
        
        curr_chs_i_sub = subsample_chs(curr_chs_i,smallest_i);
        
        subsample_regions(logical(curr_chs_i_sub)) = curr_num_i;
    end
    
    %% Show subsample vs original
    if 0
    loc_showo = cell(length(curr_regions),1);
    lat_showo = cell(length(curr_regions),1);
    loc_showo(~isnan(curr_regions)) = locs(curr_regions(~isnan(curr_regions)));
    lat_showo(~isnan(curr_regions)) = lats(curr_regions(~isnan(curr_regions)));    
             
    loc_show = cell(length(curr_regions),1);
    lat_show = cell(length(curr_regions),1);
    loc_show(~isnan(subsample_regions)) = locs(subsample_regions(~isnan(subsample_regions)));
    lat_show(~isnan(subsample_regions)) = lats(subsample_regions(~isnan(subsample_regions)));
    table(loc_showo,lat_showo,loc_show,lat_show)
    pause
    end
    
    
    % Loop over regions
    for ir = 1:nregions
        
        % get this region num
        curr_num_i = atlas_nums(ir);
        
        % get the channels in this region
        curr_chs_i = ismember(curr_regions,curr_num_i);
          
        % soz (1 if any electrodes in that region are SOZ)
        soz_atlas(ir,ip) = any(curr_soz(curr_chs_i)==1);
        
        % spikes (mean spike rate of channels in that region)
        spike_atlas(ir,ip) = nanmean(curr_spikes(curr_chs_i));
        
        for jr = 1:ir-1
            
            % get the channels in this r
            curr_num_j = atlas_nums(jr);
            curr_chs_j = ismember(curr_regions,curr_num_j);
            
          
            %% If applicable, randomly subsample electrodes to smallest number
            if force_same_num
                orig_curr_chs_i = curr_chs_i;
                % get the (subsampled) channels in this region
                curr_chs_i = ismember(subsample_regions,curr_num_i);
                curr_chs_j = ismember(subsample_regions,curr_num_j);
            end
            
            
            % get fc -  mean fc of channels in region i and channels in
            % region j
            fc_atlas(ir,jr,ip) = nanmean(curr_fc(curr_chs_i,curr_chs_j),'all');
            fc_atlas(jr,ir,ip) = nanmean(curr_fc(curr_chs_i,curr_chs_j),'all'); % symmetric matrix
            
            coh_atlas(ir,jr,:,ip) = nanmean(curr_coh(curr_chs_i,curr_chs_j,:),[1 2]);
            coh_atlas(jr,ir,:,ip) = nanmean(curr_coh(curr_chs_i,curr_chs_j,:),[1 2]);
            
        end
        
    end
    
end

soz_atlas = logical(soz_atlas);

end