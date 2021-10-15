function [soz_labels,soz_chs] = parse_soz(text,labels,name)

if isempty(text)
    soz_labels = {};
    soz_chs = [];
    return
end

% split up the labels
C = strsplit(text,',');

% clean them
soz_labels = decompose_labels(C,name);

% find the indices
[~,soz_chs] = ismember(soz_labels,labels);

if any(soz_chs == 0)
    fprintf('\n\nWARNING, cannot find the following contacts for %s:\n',name);
    soz_labels((soz_chs==0))
    
    fprintf('\nThe electrode names for this patient are:\n')
    labels
    
    fprintf('\nCheck for funny label names.\n');
end


end