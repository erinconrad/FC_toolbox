function add_file_start_times

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
script_folder = locations.script_folder;

%% Get pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Addpath
addpath(genpath(script_folder));

%% Get the xls file with the file start times
T = readtable('Manual validation.xlsx','Sheet','FileStartTimes');

%% Loop over patients
for p = 1:length(pt)
    name = pt(p).name;
    
    % find corresponding row in table
    row = strcmp(name,T.Var2);
    assert(sum(row) == 1);
    
    % Loop over files
    for f = 1:length(pt(p).ieeg.file)
        
        % get the start time of this file
        file_column_name = sprintf('Var%d',f+2); % file 1 is Var3
        start_time = T.(file_column_name)(row);
        
        assert(~isnan(start_time));
        
        pt(p).ieeg.file(f).start_time = start_time;
        
    end
    
end

%% Look for missing pts
for p = 1:length(pt)
    for f = 1:length(pt(p).ieeg.file)
        if ~isfield(pt(p).ieeg.file(f),'start_time') || isnan(pt(p).ieeg.file(f).start_time)
            error('what');
        end
        
    end
end

save([data_folder,'pt.mat'],'pt');

end