function define_run_times

%% Establish parameters
overwrite = 0;
% Decide frequency and duration of sampling
mini_block = 60; % each run should be 60 seconds
block = 600; % sample every 600 seconds

%% Seed random number generator (so that I get the same thing every time if I re-run)
rng(0)

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;


for p = 1:length(pt)

    if isempty(pt(p).ieeg)
        continue
    end
    
    if isfield(pt(p).ieeg.file(1),'run_times') == 1 && ...
            ~isempty(pt(p).ieeg.file(1).run_times)
        if overwrite == 0
            fprintf('\nskipping %s\n',pt(p).name);
        end
        continue
    end
    
    fprintf('\nDoing %s\n',pt(p).name);
    
    nfiles = length(pt(p).ieeg.file);
    for f = 1:nfiles

        % get duration
        dur = pt(p).ieeg.file(f).duration;
        
        % Split duration into chunks
        nblocks = ceil(dur/block);
        bs = (0:block:floor(dur))';
        be = bs + block;
        be(end) = dur;
        pt(p).ieeg.file(f).block_times = [bs,be];
        pt(p).ieeg.file(f).run_times = nan(nblocks,2);
        
        for b = 1:nblocks
            s = randi([0 block-mini_block]);
            start_time = min(bs(b) + s,be(b) - mini_block);
            end_time = start_time + mini_block;
            run = [start_time,end_time];
            pt(p).ieeg.file(f).run_times(b,:) = run;
        end
        
    end
    
    % Save
    save([data_folder,'pt.mat'],'pt');
end

end