function sleep_seeg_validation

% I did this for HUP217 and the general W-S order agrees with time of day
% and I looked at the scalp for a couple of random points and it generally
% aligned....so this looks promising.

% I need to figure out how to compare my manual markings with these.

% I also need to run this for more patients...

ref_start_time = datetime('01/01/2000 00:00:00','InputFormat','MM/dd/yyyy hh:mm:ss');

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_out_dir = [results_folder,'edf_out/'];
sleep_manual_dir = [results_folder,'analysis/sleep/erin_designations/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% initialize comparison array
comparison_array = {};

% Loop over folders in edf out
listing = dir([sleep_manual_dir,'*mat']);
for l = 1:length(listing)
    
    fname = listing(l).name;
    pt_name = strrep(fname,'.mat','');

   
    if ~exist([edf_out_dir,pt_name,'/sleep_stage.mat'])
        continue
    end

    fprintf('\nDoing %s\n',pt_name);

    % Load sleepseeg designation file
    sout = load([edf_out_dir,pt_name,'/sleep_stage.mat']);
    sout = sout.sout;

    % Get times of sleep transitions
    st_dates= sout.Summary(2:end,2);
    st_times = sout.Summary(2:end,3);
    stage = sout.Summary(2:end,4);
    
    new_st_dates = st_dates;
    for i = 1:length(new_st_dates) 
        if strcmp(new_st_dates{i},' ') % replace empty with last date
            new_st_dates{i} = new_st_dates{i-1};
        end
    end
    st_dates = new_st_dates;
    
    % Combine dates and times
    dts = cellfun(@(x,y) [x, ' ',y],st_dates,st_times,'uniformoutput',false);

    % convert to seconds into ieeg file
    dt = cellfun(@(x) datetime(x,'InputFormat','dd-MMM-yyyy HH:mm:ss'),dts);
    seeg_secs = seconds(dt-ref_start_time);

    if 0
        table(arrayfun(@(x) sprintf('%1.1f',x),seeg_secs,'uniformoutput',false),stage)
    end

    % Load the erin designation file
    nout = load([sleep_manual_dir,listing(l).name]); 
    if isfield (nout,'nout')
        nout = nout.nout;
    else
        nout = nout.out;
    end
    f = 1;
    blocks = nout.file(f).blocks;
    nblocks = length(blocks);
    for b = 1:nblocks
        run_times = nout.file(f).block(blocks(b)).run_times;
        erin = nout.file(f).block(blocks(b)).erin;

        % Figure out if first run time falls between two transitions
        run_start_minus_seeg = run_times(1)-seeg_secs;

        if all(run_start_minus_seeg>0) % this block is too late, stop loop
            break
        end

        if all(run_start_minus_seeg<0) % this block too early, go to next loop
            continue
        end

        % find the first seeg transition that comes BEFORE
        prior_transition = find(run_start_minus_seeg<0);
        prior_transition = prior_transition(1)-1;

        % confirm that run time end comes before the next one
        if run_times(2) > seeg_secs(prior_transition+1)
            continue
        end

        % if made it to here, then the run time falls between two state
        % transitions, and so I should be able to compare
        comparison_array = [comparison_array;...
            pt_name,erin,stage(prior_transition)];

    end
    

end


end