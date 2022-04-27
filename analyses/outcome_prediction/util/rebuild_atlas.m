function [fc_atlas,spike_atlas,soz_atlas,coh_atlas] = ...
    rebuild_atlas(fc,spikes,atlas_labels,regions,atlas_nums,labels,soz,coh)

npts = length(fc);
nregions = length(atlas_nums);
nfreqs = size(coh{1},3);

%% Remove ekg from atlas
for ip = 1:npts
    ekg = find_non_intracranial(atlas_labels{ip});
    atlas_labels{ip}(ekg) = [];
    regions{ip}(ekg) = [];
    
    assert(isequal(atlas_labels{ip},labels{ip}));
end

fc_atlas = nan(nregions,nregions,npts);
spike_atlas = nan(nregions,npts);
soz_atlas = zeros(nregions,npts);
coh_atlas = nan(nregions,nregions,nfreqs,npts);

for ip = 1:npts
    
    curr_spikes = spikes{ip};
    curr_fc = fc{ip};
    curr_regions = regions{ip};
    curr_soz = soz{ip};
    curr_coh = coh{ip};
    
    % Loop over regions
    for ir = 1:nregions
        
        % get this region num
        curr_num_i = atlas_nums(ir);
        
        % get the channels in this region
        curr_chs_i = ismember(curr_regions,curr_num_i);
        
        % soz
        soz_atlas(ir,ip) = any(curr_soz(curr_chs_i)==1);
        
        % spikes
        spike_atlas(ir,ip) = nanmean(curr_spikes(curr_chs_i));
        
        for jr = 1:ir-1
            
            % get the channels in this r
            curr_num_j = atlas_nums(jr);
            curr_chs_j = ismember(curr_regions,curr_num_j);
            
            % get fc
            fc_atlas(ir,jr,ip) = nanmean(curr_fc(curr_chs_i,curr_chs_j),'all');
            fc_atlas(jr,ir,ip) = nanmean(curr_fc(curr_chs_i,curr_chs_j),'all');
            
            coh_atlas(ir,jr,:,ip) = nanmean(curr_coh(curr_chs_i,curr_chs_j,:),[1 2]);
            coh_atlas(jr,ir,:,ip) = nanmean(curr_coh(curr_chs_i,curr_chs_j,:),[1 2]);
            
        end
        
    end
    
end

soz_atlas = logical(soz_atlas);

end