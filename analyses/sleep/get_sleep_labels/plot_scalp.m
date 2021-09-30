function plot_scalp(pt)

%% Parameters
overwrite = 0;

% sample every how many blocks
block_stride = 24; % 4 hours = 24 ten-minute blocks

% Ignore this number of blocks at beginning and end of file
buffer = 36;% 6 hours = 36 ten-minute blocks

% scalp channels I want
scalp_labels = {'F3';'C3';'FZ';'CZ';'F4';'C4';'F7';'F8'};
montage = [1 2;3 4;5 6];

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
out_dir = [results_folder,'analysis/sleep/test_data/'];
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

% ieeg stuff
login_name = locations.ieeg_login;
ieeg_folder = locations.ieeg_folder;
pwfile = locations.ieeg_pw_file;

% add script and ieeg folder to path
addpath(genpath(scripts_folder));
addpath(genpath(ieeg_folder));

%% Find pts with scalp
indices = find_pts_with_scalp(pt);

%% Loop over pts
for idx = 1:length(indices)
    p = indices(idx);
    name = pt(p).name;
    out_name = [name,'.mat'];
    
    % See if it already exists
    if overwrite == 0
        if exist([out_dir,out_name],'file') ~= 0
            out = load([out_dir,out_name]);
            out = out.out;
            
            last_file = out.last_file;
            last_run = out.last_run;
            finished = out.finished;
            
            if finished == 1
                % We have finished already
                fprintf('\nAlready finished %s, skipping...\n',name)
                continue;
            end
        else
            last_file = 0;
            last_run = buffer-block_stride;
            out.name = name;
            out.finished = 0;
            out.file(1).blocks = [];
        end
    else
        last_file = 0;
        last_run = buffer-block_stride; % ok to redo one
        pc.name = name;
        out.finished = 0;
        out.file(1).blocks = [];
    end
    
    fprintf('\nDoing %s\n',name);
     
    % loop over files
    for f = last_file+1:length(pt(p).ieeg.file)
        
        fname = pt(p).ieeg.file(f).name;
        fs = pt(p).ieeg.file(f).fs;
        
        %% Check if I have the channels I want
        chLabels = pt(p).ieeg.file(f).chLabels;
        chLabels = decompose_labels(chLabels,name);
        a = ismember(chLabels,scalp_labels);
        if sum(a) < length(scalp_labels)
            fprintf('\nmissing some scalp labels,skipping patient\n');
            skip_pt = 1;
            break
        end
        
        nblocks = size(pt(p).ieeg.file(f).block_times,1);
        
        %% loop over blocks (skipping buffer, striding)
        for b = last_run+block_stride:block_stride:nblocks-buffer
            
            fprintf('\nDoing %s file %d of %d, block %d of %d (stride %d)\n',...
                name,f,length(pt(p).ieeg.file),b,nblocks,block_stride);
            
            % Get run times
            run_times = pt(p).ieeg.file(f).run_times(b,:);
            
            %% Download ieeg data for select chs
            data = download_ieeg_select_chs(fname, login_name, pwfile, ...
            run_times,scalp_labels,name);
            values = data.values;
            
            %% Filters
            values = bandpass_filter(values,fs);
            values = notch_filter(values,fs);
            
            %% Get scalp montages
            [bi_values,bi_labels,trans_values,trans_labels] = ...
                scalp_montages(values,scalp_labels);
        
            % Remove the empty channels
            empty_chs = cellfun(@(x) strcmp(x,'-'),bi_labels);
            bi_values(:,empty_chs) = [];
            bi_labels(empty_chs) = [];
            
            empty_chs = cellfun(@(x) strcmp(x,'-'),trans_labels);
            trans_values(:,empty_chs) = [];
            trans_labels(empty_chs) = [];
            
            % Combine these channels
            all_values = [bi_values,trans_values];
            all_labels = [bi_labels;trans_labels];
            
            %% Plot
            if 0
                show_scalp_eeg(all_values,fs,all_labels)
                pause
                close all
            end
            
            %% Save the data
            out.last_file = f-1;
            out.last_run = b;
            out.name = name;
            out.file(f).block(b).values = all_values;
            out.file(f).block(b).labels = all_labels;
            out.file(f).block(b).run_times = run_times;
            out.file(f).blocks = [out.file(f).blocks;b];
            
            if b == nblocks-buffer
                out.finished = 1;
            end
            
            save([out_dir,out_name],'out');
            
        end
        
        % at end of file, reset last run
        last_run = buffer-block_stride;
    end
    
end


end