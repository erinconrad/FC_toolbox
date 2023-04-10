function T = find_stim_annotations(whichPts)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

if isempty(whichPts)
    whichPts = 1:length(pts);
end

poss_stim_text = {};
poss_stim_start = [];
files = [];
names = [];
which_ann = [];
which_event = [];

for p = whichPts
    for f = 1:length(pt(p).ieeg.file)

         if ~isfield(pt(p).ieeg.file(f),'ann'), continue; end
         n_anns = length(pt(p).ieeg.file(f).ann);
    
        
    
        for a = 1:n_anns
            
            if strcmp(pt(p).ieeg.file(f).ann,'empty'), continue; end
            
            n_events = length(pt(p).ieeg.file(f).ann(a).event);
            for i = 1:n_events
    
    
                description = pt(p).ieeg.file(f).ann(a).event(i).description;

                if contains(description,'stim','IgnoreCase',true) || ...
                        contains(description,'map','IgnoreCase',true) || ...
                        contains(description,'ms','IgnoreCase',true)
                    poss_stim_text = [poss_stim_text;description];
                    poss_stim_start = [poss_stim_start;pt(p).ieeg.file(f).ann(a).event(i).start];
                    files = [files;f];
                    names = [names;pt(p).name];
                    which_ann = [which_ann;a];
                    which_event = [which_event;i];

                end


            end

        end


    end

end

T = table(names,files,poss_stim_start,poss_stim_text,which_ann,which_event);


end