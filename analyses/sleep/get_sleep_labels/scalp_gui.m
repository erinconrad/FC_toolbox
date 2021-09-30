function scalp_gui

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
out_dir = [results_folder,'analysis/sleep/test_data/'];

% add script and ieeg folder to path
addpath(genpath(scripts_folder));

%% Loop over output files
listing = dir([out_dir,'*.mat']);
for l = 1:length(listing)
    fname = listing(l).name;
    out = load([out_dir,fname]);
    out = out.out;
    
    if isfield(out.erin_done) && out.erin_done == 1
        continue
    end
    
    for f = 1:length(out.file)
        for b = blocks
            
            if isfield(out.file(f).block(b).erin) && ~isempty(out.file(f).block(b).erin)
                values = out.file(f).block(b).values;
                labels = out.file(f).block(b).labels;
                
                
            end
            
        end
        
    end
    
end

end