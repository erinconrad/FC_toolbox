function [clean_labels,elecs,numbers] = decompose_labels(chLabels)

%{
This function takes an arbitrary set of electrode labels. It returns clean_labels,
which contains simplified labels removing leading zeros, and other bonus
parts. It also returns, for each label, the electrode name and the number
of the contact on that electrode
%}

clean_labels = cell(length(chLabels),1);
elecs = cell(length(chLabels),1);
numbers = nan(length(chLabels),1);

for ich = 1:length(chLabels)
    if ischar(chLabels)
        label = chLabels;
    else
        label = chLabels{ich};
    end

    %% if it's a string, convert it to a char
    if strcmp(class(label),'string')
        label = convertStringsToChars(label);
    end
    
    %% Remove leading zero
    % get the non numerical portion
    label_num_idx = regexp(label,'\d');
    if ~isempty(label_num_idx)

        label_non_num = label(1:label_num_idx-1);

        label_num = label(label_num_idx:end);

        % Remove leading zero
        if strcmp(label_num(1),'0')
            label_num(1) = [];
        end

        label = [label_non_num,label_num];
    end
    
    %% Remove 'EEG '
    eeg_text = 'EEG ';
    if contains(label,eeg_text)
        eeg_pos = regexp(label,eeg_text);
        label(eeg_pos:eeg_pos+length(eeg_text)-1) = [];
    end
    
    %% Remove '-Ref'
    ref_text = '-Ref';
    if contains(label,ref_text)
        ref_pos = regexp(label,ref_text);
        label(ref_pos:ref_pos+length(ref_text)-1) = [];
    end
    
    %% Remove spaces
    if contains(label,' ')
        space_pos = regexp(label,' ');
        label(space_pos) = [];
    end
    
    %% Remove '-'
    label = strrep(label,'-','');
    
    %% Fill the clean label
    clean_labels{ich} = label;
    
    %% Get the non-numerical portion
    label_num_idx = regexp(label,'\d');
    label_non_num = label(1:label_num_idx-1);
    elecs{ich} = label_non_num;
    
    %% get numerical portion
    label_num = str2num(label(label_num_idx:end));
    if isempty(label_num), label_num = nan; end
    numbers(ich) = label_num;
    
end


end