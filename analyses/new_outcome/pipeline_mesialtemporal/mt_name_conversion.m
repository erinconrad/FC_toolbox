function out_labels = mt_name_conversion(labels,name)

out_labels = labels;

%% To avoid screwing up other patients, only do this for HUP grid/strip/depth patients
if ~contains(name,'HUP'), return; end
% Get numerical portion
N = regexp(name,'\d*','Match');
N = str2num(N{1});
if N >= 127, return; end

%% Assume it's a HUP G/S/D patient. Start translating names

% Find amygdala-ish electrodes
amy = (contains(labels,'DA') | contains(labels,'AMY')) & ~contains(labels,'DAS');

% Find hippocampal-ish electrodes
hipp = contains(labels,'DH') | contains(labels,'HIP');

% Make the new A or B label
mid = cell(length(labels),1);
mid(amy==1) = {'A'};
mid(hipp==1) = {'B'};

% Get the first letter of each electrode (laterality)
lat = cellfun(@(x) x(1),labels,'UniformOutput',false);

% get numerical portion
num = cellfun(@(x) (regexp(x,'\d*','Match')),labels,'UniformOutput',false);
num(cellfun(@isempty,num)) = {{'0'}};
num = cellfun(@(x) x{1},num,'UniformOutput',false);
num = cellfun(@(x) num2str(str2num(x)),num,'uniformoutput',false);

% Put the electrode names back together
out_labels(amy|hipp) = cellfun(@(x,y,z) [x,y,z],lat(amy|hipp),mid(amy|hipp),num(amy|hipp),'uniformoutput',false);

if 0
    table(labels,out_labels)

end

end