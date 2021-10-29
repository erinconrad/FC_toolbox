function add_box_portal_times

%% Parameters
left_pad = 60*1; % assume seizure starts 60 seconds before
right_pad = 60*10; % assume seizure lasts 10 minutes after portal time says (long time given that this might be start of clip)

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];
box_times_folder = [data_folder,'sz_times/Portal sz times/'];
box_times_file = [box_times_folder,'sz_times.xlsx'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
npts = length(pt);

%% Load the sheetnames
sheets = sheetnames(box_times_file);

%% Loop over patients in the pt mat
for p = 1:npts
    
    name = pt(p).name;
    
    % Look and see if seizure times already exist (I should skip this if
    % any seizure times already exist because I think this method is the
    % least reliable)
    nfiles = length(pt(p).ieeg.file);
    already_done = 0;
    for f = 1:nfiles
        if isfield(pt(p).ieeg.file(f),'sz_time_source')
            already_done = 1;
            break
        end
    end
    if already_done, continue; end
    
    % See if it is in the excel file
    in_excel = ismember(name,sheets);
    if ~in_excel, continue; end
    
    % read appropriate table
    T = readtable(box_times_file,'Sheet',name,'ReadVariableNames',false);
    
    % Loop over seizures
    nrows = size(T,1);
    
    for s = 1:nrows
        time = T.Var3(s);
        if iscell(time)
            time = time{1};
            if strcmp(time,'N/A'), continue; end
            time = str2num(time);
        end
        

        sz_start = time - left_pad;
        sz_end = time + right_pad;
        if ismember('Var4',T.Properties.VariableNames)
            if ~contains(T.Var4{s},'D0')
                continue; 
            end
            sz_file = str2num(T.Var4{s}(end));
        else
            sz_file = 1;
        end
        
        % fill up the info
        if ~isfield(pt(p).ieeg.file(f),'sz_times')
            pt(p).ieeg.file(sz_file).sz_times = [sz_start sz_end];
        else
            pt(p).ieeg.file(sz_file).sz_times = [...
                pt(p).ieeg.file(sz_file).sz_times;sz_start sz_end];
        end
        
        pt(p).ieeg.file(sz_file).sz_time_source = 'portal times, worst source';
        
    end
   
    
end

save([data_folder,'pt.mat'],'pt');

end