function move_erin_designations_smaller_file

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
out_dir = [results_folder,'analysis/sleep/test_data/'];
new_out_dir = [results_folder,'analysis/sleep/erin_designations/'];

if ~exist(new_out_dir,'dir'), mkdir(new_out_dir); end

% add script and ieeg folder to path
addpath(genpath(scripts_folder));

%% Loop over output files
listing = dir([out_dir,'*.mat']);
for l = 1:length(listing)
    fname = listing(l).name;
    out = load([out_dir,fname]);
    out = out.out;
    
    nout.name = out.name;
    
    if ~isfield(out,'erin_done') || out.erin_done ~= 1
        continue
    end
    
    for f = 1:length(out.file)
        
        fs = out.file(f).fs;
        nout.file(f).fs = fs;
        nout.file(f).blocks = out.file(f).blocks;
        
        for ib = 1:length(out.file(f).blocks)
            b = out.file(f).blocks(ib);
            
            nout.file(f).block(b).run_times = out.file(f).block(b).run_times;
            nout.file(f).block(b).erin = out.file(f).block(b).erin;
            
            
            
        end
        
    end
    
    out = nout;
    save([new_out_dir,fname],'out');
end