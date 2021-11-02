function add_clinician_seizure_times

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Get sheet names from sz table
T = readtable([data_folder,'sz_times/single_sheet_seizure_times.xlsx']);
nrows = size(T,1);

curr_name = '';

%% Loop over rows
for r = 1:nrows
    
    % get the filename
    filename = T.PatientID{r};
    
    
    if isempty(filename), continue; end
    
    match = 0;
    
    
    %% see if the filename is in one of the patient things
    for p = 1:length(pt)
        
        if match == 1
            break
        end
         
        nfiles = length(pt(p).ieeg.file);
        for f = 1:nfiles
            currfname = pt(p).ieeg.file(f).name;
            
            % skip if no match
            if ~strcmp(filename,currfname)
                continue
            end
            
            match = 1;

            
            %% Add sz info
            if ~isfield(pt(p).ieeg.file(f),'sz_times')
                pt(p).ieeg.file(f).sz_times = [];
            end
            
            %% Overwrite conflicting sz times
            if ~strcmp(filename,curr_name)
                pt(p).ieeg.file(f).sz_times = [];
            end
            
            
            % preference for EEC
            if ~isnan(T.EEC(r))
                sz_start = T.EEC(r);
            elseif ~isnan(T.UEO(r))
                sz_start = T.UEO(r);
            else
                continue % skip if missing
            end
            
            
            if ~isnan(T.EndTimeOnPortalInSeconds(r))
                sz_end = T.EndTimeOnPortalInSeconds(r);
            else
                continue % skip if missing
            end
            
            pt(p).ieeg.file(f).sz_times = [pt(p).ieeg.file(f).sz_times;...
                sz_start sz_end];
            pt(p).ieeg.file(f).sz_time_source = T.Epileptologist(r);
            curr_name = filename;
            
            break
        end
        
        if match == 1
            break
        end
        
    end
    
end

save([data_folder,'pt.mat'],'pt');

end