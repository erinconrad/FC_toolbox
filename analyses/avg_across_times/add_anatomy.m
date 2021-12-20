function out = add_anatomy(out)


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

for i = 1:length(out)
    out_name = out(i).name;
    
    %% Load the correct summ file
    summ = load([int_folder,out_name,'.mat']);
    summ = summ.summ;
    
    
    summ_elecs = summ.labels;
    summ_anatomy = summ.anatomy;
    
    out_elecs = out(i).labels;
    
    [lia,locb] = ismember(out_elecs,summ_elecs);
    out_anatomy = cell(length(out_elecs),1);
    out_anatomy(lia) = summ_anatomy(locb(lia));
    
    assert(isequal(out_elecs(lia),summ_elecs(locb(lia))));
    out(i).anatomy = out_anatomy;
    
end

end