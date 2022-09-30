function out = convert_to_atlas(thing,atlas,atlas_names)

npts = length(thing);

unique_atlas_names = (atlas_names);
nregions = length(unique_atlas_names);

%% Figure out the dimensions of thing
thing1 = thing{1};
assert(size(thing1,1)==length(atlas{1})) % first dimension should always be nch
if ndims(thing1) == 3 % I only think this would happen for nch x nch x freq
    which_dim = 3;
    multivariate = 1;
elseif size(thing1,2) == length(atlas{1}) % nch x nch
    which_dim = 2;   
    multivariate = 1;
elseif size(thing1,2) > 1 % I think this should only happen for nch x freq
    which_dim = 2;
    multivariate = 0;
else % should only be nch x 1
    which_dim = 1;   
    multivariate = 0;
end

%% initialize out
switch which_dim
    case 1
        out = nan(npts,nregions);
    case 2
        if multivariate == 1
            out = nan(npts,nregions,nregions);
        else
            last_dim = size(thing1,2);
            out = nan(npts,nregions,last_dim);
        end
    case 3
        last_dim = size(thing1,3);
        out = nan(npts,nregions,nregions,last_dim);
end

%% Fill it up
% Loop over pts
for ip = 1:npts
    
    curr_thing = thing{ip};
    curr_atlas = atlas{ip};

    % Loop over regions
    for ir = 1:nregions
        
        % get current atlas region name
        curr_region = unique_atlas_names{ir};

        % find matching elecs
        matching_elecs = strcmp(curr_atlas,curr_region);

        % Take average of stuff
        switch which_dim
            case 1
                out(ip,ir) = nanmean(curr_thing(matching_elecs));
            otherwise
                if ~multivariate
                    out(ip,ir,:) = nanmean(curr_thing(matching_elecs,:),1);
                else
                    for jr = 1:ir-1
                        curr_region_j = unique_atlas_names{jr};
                        matching_elecs_j = strcmp(curr_atlas,curr_region_j);
                        
                        switch which_dim
                            case 2
                                out(ip,ir,jr) = nanmean(curr_thing(matching_elecs,matching_elecs_j),'all');
                                out(ip,jr,ir) = nanmean(curr_thing(matching_elecs,matching_elecs_j),'all');
                            case 3
                                out(ip,ir,jr,:) = nanmean(curr_thing(matching_elecs,matching_elecs_j,:),[1 2]);
                                out(ip,jr,ir,:) = nanmean(curr_thing(matching_elecs,matching_elecs_j,:),[1 2]);
                        end
    
                    end
                end
            
        end

    end
end



end