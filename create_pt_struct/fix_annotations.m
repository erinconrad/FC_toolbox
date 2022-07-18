function fix_annotations

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
    pt(p).ieeg.file = rmfield(pt(p).ieeg.file,'ann'); 
    for f = 1:length(pt(p).ieeg.file)
        ieeg_name = pt(p).ieeg.file(f).name;
        session = IEEGSession(ieeg_name,login_name,pwfile);
        
        clear ann
        % Add annotations
        n_layers = length(session.data.annLayer);

        if n_layers == 0
            pt(p).ieeg.file(f).ann = 'empty';
            continue
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
            pt(p).ieeg.file(f).ann(ai) = ann;
        end
        
    end
    
end

save([data_folder,'pt.mat'],'pt');

end