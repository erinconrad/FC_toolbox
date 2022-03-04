function [SSE,idx] = find_optimal_k(thing,ks_to_check)

nk = length(ks_to_check);
nattempts = 30;


SSE = zeros(nk,1);
idx = cell(nk,1);

for ik = 1:nk
    k = ks_to_check(ik);
    SSE_temp = zeros(nattempts,1);
    idx_test_all = cell(nattempts,1);
    
    for j = 1:nattempts
        
        [idx_test,C_test] = kmeans(thing,k);
        idx_test_all{j} = idx_test;
        
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
    [SSE(ik),I] = min(SSE_temp);
    idx{ik} = idx_test_all{I};
    
end


end