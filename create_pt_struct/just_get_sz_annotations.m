function T = just_get_sz_annotations(hup_nos)

%% Get file locs
locations = fc_toolbox_locs;
ieeg_folder = locations.ieeg_folder;
script_folder = locations.script_folder;
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Add paths
addpath(genpath(ieeg_folder));
addpath(genpath(script_folder));

p = 1;

for i = 1:length(hup_nos)

    % Say I have not found the patient
    found_pt = 0;

    curr_hup_no = hup_nos(i);
    base_ieeg_name = sprintf('HUP%d_phaseII',curr_hup_no);

    % Initialize stuff for ieeg files
    dcount = 0; % which file
    add_it = 0; % should I add info
    finished = 0; % done with pt
    
    while 1
        if dcount == 0
            
            % Try to get ieeg file with just the base name
            ieeg_name = base_ieeg_name;
            attempt = 0;
            while 1
                try
                    session = IEEGSession(ieeg_name,login_name,pwfile);
                    finished = 1;
                    add_it = 1;
                    dcount = 1;
                    fprintf('\nFound %s, getting info\n',ieeg_name);
                    break
                    
    
                catch ME
                    if contains(ME.message,'503') || contains(ME.message,'504') || ...
                            contains(ME.message,'502') || contains(ME.message,'500')
                        attempt = attempt + 1;
                        fprintf('Failed to retrieve ieeg.org data, trying again (attempt %d)\n',attempt); 
                    else
                        ME
                        fprintf('\nDid not find %s, adding an appendage\n',ieeg_name);
                        if exist('session','var') ~= 0
                            session.delete;
                        end                    
                    end
            
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
            
            clear ann
            for ai = 1:n_layers
                
                count = 0;
                
                while 1 % while loop because ieeg only lets you pull 250 at once
                    clear event % important!
                    
                    % ask it to pull next (up to 250) events after count
                    if count == 0
                        a=session.data.annLayer(ai).getEvents(count);
                    else
                        a=session.data.annLayer(ai).getNextEvents(a(n_ann));
                    end
                    if isempty(a), break; end
                    n_ann = length(a);
                    for k = 1:n_ann
                        event(k).start = a(k).start/(1e6);
                        event(k).stop = a(k).stop/(1e6); % convert from microseconds
                        event(k).type = a(k).type;
                        event(k).description = a(k).description;
                    end
                    ann.event(count+1:count+k) = event;
                    
                    count = count + n_ann;
                end
                
                ann.name = session.data.annLayer(ai).name;
                pt(p).ieeg.file(dcount).ann(ai) = ann;
                
            end
            found_pt = 1; % say I found the patient
            pt(p).name = sprintf('HUP%d',curr_hup_no);
            
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
    end
            
    % advance pt index if I did find it
    if found_pt == 1
        p = p + 1;
    end
    
           
end

%% Mine annotations for szs
poss_sz_text = {};
poss_sz_start = [];
files = [];
names = [];
which_ann = [];
which_event = [];
for p = 1:length(pt)
for f = 1:length(pt(p).ieeg.file)
    
    if ~isfield(pt(p).ieeg.file(f),'ann'), continue; end
    n_anns = length(pt(p).ieeg.file(f).ann);

    

    for a = 1:n_anns
        
        if strcmp(pt(p).ieeg.file(f).ann,'empty'), continue; end
        
        n_events = length(pt(p).ieeg.file(f).ann(a).event);
        for i = 1:n_events


            description = pt(p).ieeg.file(f).ann(a).event(i).description;

            % search for seizure-y strings
            if contains(description,'seizure','IgnoreCase',true) || ...
                    contains(description,'sz','IgnoreCase',true) || ...
                    contains(description,'onset','IgnoreCase',true) || ...
                    contains(description,'UEO','IgnoreCase',true) || ...
                    contains(description,'EEC','IgnoreCase',true) 

                poss_sz_text = [poss_sz_text;description];
                poss_sz_start = [poss_sz_start;pt(p).ieeg.file(f).ann(a).event(i).start];
                files = [files;f];
                names = [names;pt(p).name];
                which_ann = [which_ann;a];
                which_event = [which_event;i];

            end
        end

    end

   
    
end
end

T = table(names,files,poss_sz_start,poss_sz_text,which_ann,which_event);

end