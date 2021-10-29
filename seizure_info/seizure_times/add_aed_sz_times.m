function add_aed_sz_times

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load Nina's AED seizure time structure
load([data_folder,'sz_times/aed_sz_times/AED_seizure_times.mat']);

nrows = size(cohort_info,1);

% Loop over rows
for r = 1:nrows
    name = cohort_info.PtIDs{r};
    file_name = cohort_info.file_used{r};
    
    
    found_pt = 0;
    
    % find the corresponding patient and file
    for p = 1:length(pt)
        if strcmp(pt(p).name,name)
            found_pt = 1;
            for f = 1:length(pt(p).ieeg.file)
                if strcmp(pt(p).ieeg.file(f).name,file_name)
                    % if match, fill up seizure info
                    sz_times = all_seizure_times{r};
                    pt(p).ieeg.file(f).sz_times = sz_times*3600; % turn from hours to seconds
                    pt(p).ieeg.file(f).sz_time_source = 'Nina AED sz times';
                    break
                end
                
            end
            
        end
        if found_pt, break; end
        if p == length(pt)
            error('did not find');
        end
    end
end

save([data_folder,'pt.mat'],'pt');

end