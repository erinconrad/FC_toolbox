function ai = get_ai(thing,last_dim,uni,elec_lats)

dims = ndims(thing{1});
npts = length(thing);

%% Get average left and right
% Initialize hemi_lr according to dimension
hemi_lr = nan(npts,2,last_dim);

% Fill it according to dimension
for ip = 1:npts
    L = strcmp(elec_lats{ip},'L'); % get electrodes on L
    R = strcmp(elec_lats{ip},'R'); % electrodes on R

    % Skip if fewer than 5 contacts on either side
    if sum(L) < 3 || sum(R) < 3
        hemi_lr(ip,:,:) = repmat([nan nan],1,1,last_dim);
        continue
    end

    switch dims
        case 2
            if isnan(uni) % elec nums
                hemi_lr(ip,:,:) = [sum(thing{ip}(L)) sum(thing{ip}(R))];
                if any(hemi_lr(ip,:)==0)
                    hemi_lr(ip,:) = [nan nan];
                end
            elseif uni && last_dim == 1 % spikes
                hemi_lr(ip,:,:) = [nanmean(thing{ip}(L)) nanmean(thing{ip}(R))];
            elseif uni && last_dim ~= 1 % bandpower
                hemi_lr(ip,1,:) = nanmean(thing{ip}(L,:),1);
                hemi_lr(ip,2,:) = nanmean(thing{ip}(R,:),1);
            elseif uni == 0 % FC (pearson)
                hemi_lr(ip,1,:) = nanmean(thing{ip}(L,L),'all');
                hemi_lr(ip,2,:) = nanmean(thing{ip}(R,R),'all');
            else
                error('what')
            end
        case 3 % coherence only?
            hemi_lr(ip,1,:) = nanmean(thing{ip}(L,L,:),[1 2]);
            hemi_lr(ip,2,:) = nanmean(thing{ip}(R,R,:),[1 2]);
    end

end

%% Calculate asymmetry index of L-R
ai = squeeze(asymmetry_index(hemi_lr(:,1,:),hemi_lr(:,2,:)));



end