function all_aucs = sleep_duration(durations,just_gray)

%% Locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

time_roc_folder = [out_folder,'time_roc/'];
if ~exist(time_roc_folder,'dir')
    mkdir(time_roc_folder)
end

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;
npts = length(out.circ_out.names);
all_soz = out.bin_out.all_is_soz;

%% Initialize stuff
ndurs = length(durations);
all_aucs = nan(2,ndurs,npts); % wake/sleep, then durations, then nbs

% Loop over wake and sleep
for iws = 1:2
    
    % Loop over durations
    for id = 1:ndurs
        fprintf('\nDoing iws %d of %d duration %d of %d\n',iws,2,id,ndurs);
        duration = durations{id};
        
        % Loop over pts
        for ip = 1:npts
            %fprintf('\npatient = %d of %d\n',ip,npts);
            curr_soz = all_soz{ip};
            if sum(curr_soz) == 0
                continue;
            else
                mout = updated_classifier_may2022(ip,1,iws,duration,just_gray);
                all_aucs(iws,id,ip) = mout.AUC;
            end
            
            
        end
        
    end
    
end

end