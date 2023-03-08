function find_intracranial_pts(ieeg_nums,site)

%% Parameters
overwrite = 0;

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
ieeg_folder = locations.ieeg_folder;
script_folder = locations.script_folder;
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Add paths
addpath(genpath(ieeg_folder));
addpath(genpath(script_folder));

%% Load data file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Remove dangling patients
npts = length(pt);
for i = npts:-1:1
    if isempty(pt(i).ieeg)
        pt(i) = [];
    else
        break
    end
end

if overwrite == 1
    p = 1;
else
    p = length(pt) + 1; % start at next index of pt struct
end

if isempty(ieeg_nums)
    first_num = 100;
    last_num = 300;
else
    first_num = ieeg_nums(1);
    last_num = ieeg_nums(2);
end

% Loop over numbers
i = first_num; % the HUP*** I am checking
while 1
    
    % Say I have not found the patient
    found_pt = 0;
    
    % Base name
    switch site
        case 'hup'
            if ismember(i,[69,71,76,77,93,94,96,97,98])
                base_ieeg_name = sprintf('HUP0%d_phaseII',i);
            else
                base_ieeg_name = sprintf('HUP%d_phaseII',i);
            end
        case 'musc'
            if i < 10
                base_ieeg_name = sprintf('MP000%d',i);
            else
                base_ieeg_name = sprintf('MP00%d',i);
            end
        case 'musc_3T'
            base_ieeg_name = sprintf('3T_MP00%d',i);
    end
    
    % Initialize stuff for ieeg files
    dcount = 0; % which file
    add_it = 0; % should I add info
    finished = 0; % done with pt
    
    while 1
        if dcount == 0
            
            % Try to get ieeg file with just the base name
            ieeg_name = base_ieeg_name;
            
            try
                session = IEEGSession(ieeg_name,login_name,pwfile);
                finished = 1;
                add_it = 1;
                dcount = 1;
            catch
                
                fprintf('\nDid not find %s, adding an appendage\n',ieeg_name);
                if exist('session','var') ~= 0
                    session.delete;
                end
                
            end
            
        else % if dcount > 0, trying appendage
            % Try it with an appendage
            ieeg_name = [base_ieeg_name,'_D0',sprintf('%d',dcount)];
            try
                session = IEEGSession(ieeg_name,login_name,pwfile);
                finished = 0;
                add_it = 1;
            catch
                add_it = 0;
                finished = 1; % if I can't find it adding appendage, nothing else to check
            end
            
        end
        
        % Add session info
        if add_it == 1
            pt(p).ieeg.file(dcount).fs = session.data.sampleRate;
            pt(p).ieeg.file(dcount).name = session.data.snapName;
            pt(p).ieeg.file(dcount).chLabels = session.data.channelLabels(:,1);
            pt(p).ieeg.file(dcount).duration = session.data.rawChannels(1).get_tsdetails.getDuration/(1e6); % convert from microseconds
            
            clear ann
            % Add annotations
            n_layers = length(session.data.annLayer);
    
            if n_layers == 0
                pt(p).ieeg.file(dcount).ann = 'empty';
            end

            for ai = 1:n_layers
                clear event
                a=session.data.annLayer(ai).getEvents(0);
                n_ann = length(a);
                for k = 1:n_ann
                    event(k).start = a(k).start/(1e6);
                    event(k).stop = a(k).stop/(1e6); % convert from microseconds
                    event(k).type = a(k).type;
                    event(k).description = a(k).description;
                end
                ann.event = event;
                ann.name = session.data.annLayer(ai).name;
                pt(p).ieeg.file(dcount).ann(ai) = ann;
            end
            
            found_pt = 1; % say I found the patient
        end
        
        % done with patient
        if finished == 1
            if exist('session','var') ~= 0
                session.delete;
            end
            break % break out of ieeg loop for that patient
        end
        
        dcount = dcount + 1; % if not finished, see if another appendage
        if exist('session','var') ~= 0
            session.delete;
        end
        
        
    end

    % done with that patient 
    switch site
        case 'hup'
            if i <100
                pt(p).name = sprintf('HUP0%d',i);
            else
                pt(p).name = sprintf('HUP%d',i);
            end
        case 'musc'
            pt(p).name = base_ieeg_name;
        case 'musc_3T'
            pt(p).name = strrep(base_ieeg_name,'3T_','');
    end
    
    % Save the file
    save([data_folder,'pt.mat'],'pt');
    
    % advance pt index if I did find it
    if found_pt == 1
        p = p + 1;
    end
    
    % advance count regardless
    i = i + 1;
    if i > last_num
        break
    end
end
