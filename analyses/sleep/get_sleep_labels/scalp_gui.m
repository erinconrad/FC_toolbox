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
    
    if isfield(out,'erin_done') && out.erin_done == 1
        continue
    end
    
    for f = 1:length(out.file)
        fs = out.file(f).fs;
        for ib = 1:length(out.file(f).blocks)
            b = out.file(f).blocks(ib);
            
            if isfield(out.file(f).block(b),'erin') && ~isempty(out.file(f).block(b).erin)
                continue
            end
            
            values = out.file(f).block(b).values;
            labels = out.file(f).block(b).labels;
            show_scalp_eeg(values,fs,labels)

            while 1
                prompt = sprintf(['\nBlock %d of %d Name the state:\n'...
                    'a = awake\n1 = N1\n2 = N2\n3 = N3\nr = REM\n'...
                    'u = unclear\nz = seizure\n'],ib,length(out.file(f).blocks));
                str = input(prompt,'s');
                if ~ismember(str,{'a','1','2','3','r','u','z'})
                    fprintf('\nRe-enter your choice.\n');
                    continue
                else
                    out.file(f).block(b).erin = str;
                    break
                end
            end
            close all
                    
            save([out_dir,fname],'out');
        end
        
    end
    out.erin_done = 1;
    
end

end