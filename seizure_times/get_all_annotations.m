function get_all_annotations

%% Parameters

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

for p = 1:length(pt)
    for f = 1:length(pt(p).ieeg.file)
        ieeg_name = pt(p).ieeg.file(f).name;
        session = IEEGSession(ieeg_name,login_name,pwfile);
        
        n_layers = length(session.data.annLayer);
        
        if n_layers == 0
            pt(p).ieeg.file(f).ann = 'empty';
        end
        
        
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
            pt(p).ieeg.file(f).ann(ai) = ann;
            
        end
        
        
    end
end

save([data_folder,'pt.mat'],'pt');

end