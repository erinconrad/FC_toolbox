function mine_annotations_for_szs

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
out_folder = [data_folder,'sz_times/'];

%% Load data file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;


poss_sz_text = {};
poss_sz_start = [];
files = [];
names = [];
for p = 1:14
for f = 1:length(pt(p).filename)
    
    if ~isfield(pt(p).filename(f),'ann'), continue; end
    n_anns = length(pt(p).filename(f).ann);

    

    for a = 1:n_anns
        
        n_events = length(pt(p).filename(f).ann(a).event);
        for i = 1:n_events


            description = pt(p).filename(f).ann(a).event(i).description;

            % search for seizure-y strings
            if contains(description,'seizure','IgnoreCase',true) || ...
                    contains(description,'sz','IgnoreCase',true) || ...
                    contains(description,'onset','IgnoreCase',true) || ...
                    contains(description,'UEO','IgnoreCase',true) || ...
                    contains(description,'EEC','IgnoreCase',true) 

                poss_sz_text = [poss_sz_text;description];
                poss_sz_start = [poss_sz_start;pt(p).filename(f).ann(a).event(i).start];
                files = [files;f];
                names = [names;pt(p).name];

            end
        end

    end

   
    
end
end

T = table(names,files,poss_sz_start,poss_sz_text);
writetable(T,[out_folder,'possible_sz_annotations.csv'])



end