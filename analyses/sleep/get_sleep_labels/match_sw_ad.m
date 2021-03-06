function summ = match_sw_ad(pt,summ)

for l = 1:length(summ)
    name = summ(l).name;
    
    % find corresponding pt p
    p = find_matching_pt(name,pt);
    if isnan(p), error('did not find pt'); end
    
    % Loop over sleep and wake tables
    fn = fieldnames(pt(p).sw);
    for n = 1:length(fn)
        T = pt(p).sw.(fn{n});
        
        % Prep array of run indices
        run_indices = nan(size(T,1),1);
        
        % Loop over annotations in the table
        for t = 1:size(T,1)
            % get file and time
            f = T{t,1};
            time = T{t,2};
            
            % find matching indices in pt
            matching_files = summ(l).file_index == f;
            
            % find the index of the closest time
            temp_file_times = summ(l).file_times;
            temp_file_times(~matching_files) = inf; % this is just to make sure I will only find closest time in matching file
            [~,I] = min(abs(temp_file_times - time));
            
            % add this to run_indices array
            run_indices(t) = I;
            
        end
        
        % Remove duplicates in run indices array
        run_indices = unique(run_indices);
        
        % get the alpha delta ratios of these run indices
        ad = summ(l).ad(:,run_indices);
        
        % average across electrodes
        ad = nanmean(ad,1);
        
        % add data to summ
        summ(l).sw.(fn{n}) = ad;
    end
    
end

end