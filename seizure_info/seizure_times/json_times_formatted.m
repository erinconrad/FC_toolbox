function T = json_times_formatted

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];

%% Load json file into struct
json_file{1} = [data_folder,'sz_times/DATA_MASTER.json'];
json_file{2} = [data_folder,'sz_times/akash_sz_times/DATA_MASTER.json'];
json_file{3} = [data_folder,'sz_times/nina_sz_times/DATA_MASTER.json'];

T = cell2table(cell(0,3));
T.Properties.VariableNames{'Var1'} = 'Name';
T.Properties.VariableNames{'Var2'} = 'Start';
T.Properties.VariableNames{'Var3'} = 'End';

% Loop over which json files to add
for js = 1:length(json_file)
    curr_json_file = json_file{js};

    data = loadjson(curr_json_file);
    ptnames = fieldnames(data.PATIENTS);
    
    for p = 1:length(ptnames)
        
        curr = data.PATIENTS.(ptnames{p});
        
        % Add seizure times
        if isfield(curr,'Events')
            if isfield(curr.Events,'Ictal')
                szs = curr.Events.Ictal;
                szfnames = fieldnames(szs);

                for s = 1:length(szfnames)
                    curr_sz = szs.(szfnames{s});
                    if isfield(curr_sz,'SeizureEEC') && ~isempty(curr_sz.SeizureEEC)
                        sz_start = curr_sz.SeizureEEC;
                    else
                        sz_start = curr_sz.SeizureUEO;
                    end
                    
                    if ~isfield(curr_sz,'SeizureEnd')
                        sz_end = nan;
                    else
                        sz_end = curr_sz.SeizureEnd;
                    end
                    
                    T = [T;{ptnames{p},sz_start,sz_end}];

                end
            end
        end
        
        
        
    end
    
end

end