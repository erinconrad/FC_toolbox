function SSE = find_optimal_k(thing)

max_k = 10;
nattempts = 30;


SSE = zeros(max_k,1);

for k = 1:max_k
    SSE_temp = zeros(nattempts,1);
    
    for j = 1:nattempts
        
        [idx_test,C_test] = kmeans(thing,k);
        
        % Get SSE
        for i = 1:k
            
            % Get those assigned to cluster
            curr_idx = find(idx_test == i);
            curr_C = C_test(i,:);
            
            % repmat C
            curr_C_rep = repmat(curr_C,length(curr_idx),1);
            
            % Add SS distance between each observation and its cluster
            % centroid
            SSE_temp(j) = SSE_temp(j) + ...
                sum((curr_C_rep - thing(curr_idx,:)).^2,'all');
            
        end
        
    end
    
    % select the best
    SSE(k) = min(SSE_temp);
    
end


end