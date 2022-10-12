function add_sztimes_from_manual_valdation_file(overwrite)

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load AllSeizureTimes sheet of Manual validation file
T = readtable('Manual validation.xlsx','Sheet','AllSeizureTimes');
nrows = size(T,1);

curr_szs = [];
curr_source = {};
curr_semiology = {};
curr_pt = 1;
curr_f = 1;

for r = 1:nrows
    
    % Get patient name
    name = T.Patient{r};
    
    if isempty(name), continue; end % I have nan rows between patients
    
    % find corresponding pt
    found_it = 0;
    for ip = 1:length(pt)
        if strcmp(pt(ip).name,name)
            found_it = 1;
            break
        end
            
    end
    
    assert(found_it)
    
    % get file
    f = T.IEEGID(r);
    
    % see if already data, skip if overwrite == 0
    if overwrite == 0
        if isfield(pt(ip).ieeg.file(f),'sz_times') && ~isempty(pt(ip).ieeg.file(f).sz_times)
            fprintf('\nSkipping %s file %d\n',name,f);
            continue;
        end
    end
    
    fprintf('\nDoing %s file %d\n',name,f);
    
    % if file or patient is new from prior, add to struct and reset
    if ip ~= curr_pt || curr_f ~= f
        
        pt(curr_pt).ieeg.file(curr_f).sz_times = curr_szs;
        pt(curr_pt).ieeg.file(curr_f).sz_time_source = curr_source;
        pt(curr_pt).ieeg.file(curr_f).sz_semiology = curr_semiology;
        
        curr_pt = ip; curr_f = f;
        curr_szs = [];
        curr_source = {};
        curr_semiology = {};
    end
       
    % get sz data
    sz_start = T.start(r);
    sz_end = T.xEnd(r);
    sz_source = T.source{r};
    sz_semiology = T.Semiology{r};
    
    % add to curr table
    curr_szs = [curr_szs;sz_start sz_end];
    curr_source = [curr_source;sz_source];
    curr_semiology = [curr_semiology;sz_semiology];
    
end

% if last row, add for last one
if ~isempty(curr_szs)
    if isfield(pt(ip).ieeg.file(f),'sz_times') && ~isempty(pt(ip).ieeg.file(f).sz_times)
            fprintf('\nSkipping %s file %d\n',name,f);
    else
        fprintf('\nDoing %s file %d\n',name,f);
    end
    pt(curr_pt).ieeg.file(curr_f).sz_times = curr_szs;
    pt(curr_pt).ieeg.file(curr_f).sz_time_source = curr_source;
    pt(curr_pt).ieeg.file(curr_f).sz_semiology = curr_semiology;
end

save([data_folder,'pt.mat'],'pt');
    
end