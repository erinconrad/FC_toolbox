function scalp_gui

%{
This is the function to generate "true" sleep/wake labels to test the alpha
delta ratio performance. It takes a bunch of minutes selected every four
hours and it plots the scalp eeg in a bipolar montage (a limited set of
electrodes). Erin then goes through and designates the state and it saves
the output in the erin_designations folder. 

Once these are made, I can then run compare_sw_ad, which takes the alpha
delta ratio in the corresponding minutes and then I build an ROC curve to
see how well the alpha delta ratio does at discriminating sleep vs wake
state

To run this, I need to move patient test_data files over to my laptop one
at a time (because they're kind of big).
%}


%% Parameters
which_montage = 'bipolar';
skip_chs = 'FZ-CZ';

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
    
    % load the test data file
    out = load([out_dir,fname]);
    out = out.out;
    name = out.name;
    
    % load the final output file
    if exist([new_out_dir,name,'.mat'],'file') ~= 0
        nout = load([new_out_dir,name,'.mat']);
        if isfield(nout,'out')
            nout = nout.out;
        elseif isfield(nout,'nout')
            nout = nout.nout;
        end
    else
        nout.name = out.name;
    end
    
    
    if isfield(out,'erin_done') && out.erin_done == 1
        continue
    end
    
    
    
    for f = 1:length(out.file)
        fs = out.file(f).fs;
        
        nout.file(f).fs = fs;
        nout.file(f).blocks = out.file(f).blocks;
        
        for ib = 1:length(out.file(f).blocks)
            b = out.file(f).blocks(ib);
            
            if isfield(nout.file(f),'block') && ...
                    length(nout.file(f).block) >=b &&...
                    isfield(nout.file(f).block(b),'erin') &&...
                    ~isempty(nout.file(f).block(b).erin)
                continue
            end
            
            if strcmp(which_montage,'bipolar')
                values = out.file(f).block(b).bi_values;
                labels = out.file(f).block(b).bi_labels;
            elseif strcmp(which_montage,'transverse')
                values = out.file(f).block(b).trans_values;
                labels = out.file(f).block(b).trans_labels;
            end
            to_skip = ismember(labels,skip_chs);
            values(:,to_skip) = [];
            labels(to_skip) = [];
            show_scalp_eeg(values,fs,labels)

            while 1
                prompt = sprintf(['\n%s file %d of %d Block %d of %d\nName the state:\n'...
                    'a = awake\n1 = N1\n2 = N2\n3 = N3\nr = REM\n'...
                    'u = unclear\nz = seizure\n'],...
                    out.name,f,length(out.file),ib,length(out.file(f).blocks));
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
            
            nout.file(f).block(b).run_times = out.file(f).block(b).run_times;
            nout.file(f).block(b).erin = out.file(f).block(b).erin;
                    
            save([new_out_dir,fname],'nout');
        end
        
    end
    out.erin_done = 1;
    save([out_dir,fname],'out');
    
end

end