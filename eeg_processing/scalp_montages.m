function [bi_values,bi_labels,trans_values,trans_labels] = ...
    scalp_montages(values,labels)

%% Define montages
bi_pair = {'Fp1','F7';...
    'F7','T3';...
    'T3','T5';...
    'T5','O1';...
    'Fp2','F8';...
    'F8','T4';...
    'T4','T6';...
    'T6','O2';...
    'Fp1','F3';...
    'F3','C3';...
    'C3','P3';...
    'P3','O1';...
    'Fp2','F4';...
    'F4','C4';...
    'C4','P4';...
    'P4','O2';...
    'FZ','CZ'};
trans_pair = {'F7','F3';...
    'F3','FZ';...
    'FZ','F4';...
    'F4','F8';...
    'T3','C3';...
    'C3','CZ';...
    'CZ','C4';...
    'C4','T4'};

%% initialize stuff
ntimes = size(values,1);
bi_values = nan(ntimes,size(bi_pair,1));
bi_labels = cell(size(bi_pair,1),1);
trans_values = nan(ntimes,size(trans_pair,1));
trans_labels = cell(size(trans_pair,1),1);

%% Build bipolar montage
for i = 1:size(bi_pair,1)
    first = bi_pair{i,1};
    second = bi_pair{i,2};
    
    % find the corresponding elements
    first_idx = (strcmp(labels,first));
    second_idx = (strcmp(labels,second));
    
    if sum(first_idx) == 0 || sum(second_idx) == 0
        bi_labels{i} = '-';
        continue
    end
    
    % define the value to be first minus second
    bi_values(:,i) = values(:,first_idx) - values(:,second_idx);
    
    % define the label
    bi_labels{i} = [first,'-',second];
    
end

%% Build trans montage
for i = 1:size(trans_pair,1)
    first = trans_pair{i,1};
    second = trans_pair{i,2};
    
    % find the corresponding elements
    first_idx = (strcmp(labels,first));
    second_idx = (strcmp(labels,second));
    
    if sum(first_idx) == 0 || sum(second_idx) == 0
        trans_labels{i} = '-';
        continue
    end
    
    % define the value to be first minus second
    trans_values(:,i) = values(:,first_idx) - values(:,second_idx);
    
    % define the label
    trans_labels{i} = [first,'-',second];
    
end


end