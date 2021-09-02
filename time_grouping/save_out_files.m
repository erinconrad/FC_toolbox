function save_out_files

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
file_folder = [results_folder,'all_out/'];
out_folder = [file_folder,'summary_files/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

listing = dir([file_folder,'*.mat']);
for l = 1:length(listing)
    fname = [file_folder,listing(l).name];
    pc = load(fname);
    pc = pc.pc;
    
    out = net_over_time(pc);
    out = reconcile_files(out);
    out.name = pc.name;
    
    out = rmfield(out,'file');
    
    save([out_folder,out.name],'out')
end




end