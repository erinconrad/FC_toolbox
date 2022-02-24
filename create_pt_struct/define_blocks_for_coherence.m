function define_blocks_for_coherence

%% Establish parameters
overwrite = 0;
block_span = 24; % span 24 10-minute blocks at a time (4 hours at a time)

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
    
    if isfield(pt(p).ieeg.file(1),'coherence_blocks') == 1 &&...
       ~isempty(pt(p).ieeg.file(1).coherence_blocks)
        if overwrite == 0
            fprintf('\nskipping %s\n',pt(p).name);
        end
        continue
    end

    fprintf('\nDoing %s\n',pt(p).name);
    
    nfiles = length(pt(p).ieeg.file);
    for f = 1:nfiles
        
        % get number of blocks
        nblocks = size(pt(p).ieeg.file(f).block_times,1);
        
        % initialize binary array nblocks x 1 stating if doing coherence
        coherence_blocks = zeros(nblocks,1);
        
        % Split into chunks
        nchunks = ceil(nblocks/block_span);
        chunk_start = (1:block_span:nblocks)';
        chunk_end =chunk_start+block_span;
        chunk_end(end) = nblocks;
        
        % Loop over chunks
        for ic = 1:nchunks
            % pick a random integer in the chunk
            rand_block = randi([chunk_start(ic) chunk_end(ic)]);
            
            % Make that block be a 1 for coherence blocks
            coherence_blocks(rand_block) = 1;
        end
        
        pt(p).ieeg.file(f).coherence_blocks = coherence_blocks;

    end
    
    % Save
    save([data_folder,'pt.mat'],'pt');

end