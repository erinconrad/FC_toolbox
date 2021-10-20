function adj_out = wrap_or_unwrap_adjacency_fc_toolbox_old(adj_in)


%% Detect if it is currently flat or expanded
[m,n] = size(adj_in);
smallest_dim = min([m,n]);
if m == n
    currently_flat = 0;
else
    currently_flat = 1;
    if smallest_dim == 1
        single_dim = 1;
        ntimes = 1;
    else
        single_dim = 0;
        ntimes = n;
    end
end


%% Expand it if it's currently flat
if currently_flat == 1
    
    % Get the number of channels (solving quadratic equation)
    
    y = length(adj_in(:,1));
    
    nchs = 0.5 + sqrt(0.25+2*y);
    
    % sanity checks
    if (nchs-floor(nchs)) > 0.01
        error('what\n');
    end
    nchs = round(nchs);
    if nchs*(nchs-1)/2 ~= y
        error('what\n');
    end
    
    adj_out = nan(nchs,nchs,ntimes);
    
    for t = 1:ntimes
        curr_adj = adj_in(:,t);

        % Reconstruct upper triangular matrix
        curr_out = zeros(nchs,nchs);
        count = 0;
        for i = 1:nchs
            for j = 1:i-1
                count = count + 1;
                curr_out(j,i) = curr_adj(count);
            end
        end

        if count ~= length(curr_adj)
            error('what\n');
        end

        % Reflect across the diagonal to get full adjacency matrix
        curr_out = curr_out + curr_out';
        
        adj_out(:,:,t) = curr_out;
    end
    
else
    
    %% Flatten it if it's currently expanded
    nchs = size(adj_in,1);
    adj_out = zeros(nchs*(nchs-1)/2,1);
    
    count = 0;
    
    for i = 1:nchs
        for j = 1:i-1
            count = count + 1;
            adj_out(count,:) = adj_in(j,i);
        end
    end
    
end

a = size(adj_out);
if length(a) == 3
    if a(3) == 1
        adj_out = squeeze(adj_out);
    end
end

end