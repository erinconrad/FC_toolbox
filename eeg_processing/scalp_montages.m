function [bi_values,bi_labels,trans_values,trans_labels] = ...
    scalp_montages(values,labels)

%% Define montages
bi_pair = {'F3','C3';...
    'FZ','CZ';
    'F4','C4'};
trans_pair = {'F7','F3';...
    'F3','FZ';...
    'FZ','F4';...
    'F4','F8';...
    'C3','CZ';...
    'CZ','C4'};

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
    
    % define the value to be first minus second
    trans_values(:,i) = values(:,first_idx) - values(:,second_idx);
    
    % define the label
    trans_labels{i} = [first,'-',second];
    
end


end