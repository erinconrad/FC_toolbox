function erase_all_sleep_labels(name)

prompt = sprintf('\nAre you sure? (y/n)\n');
str = input(prompt,'s');
if ~strcmp(str,'y')
    return
end

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
out_dir = [results_folder,'analysis/sleep/test_data/'];


out = load([out_dir,name,'.mat']);
out = out.out;

for f = 1:length(out.file)
        for ib = 1:length(out.file(f).blocks)
            b = out.file(f).blocks(ib);
            out.file(f).block(b).erin = [];
            
        end
        
end

save([out_dir,name,'.mat'],'out')

end