function auc_duration(nb)

%% Parameters
%durations = {30,60*4,[]};
durations = {1, 5, 10, 30, 60, 60*2,[]};
myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];

%% Seed rng (for splitting testing and training data)
% If I don't seed this the AUC usually bounces around 0.73-0.82
rng(0)

%% Locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
time_roc_folder = [out_folder,'time_roc/'];
if ~exist(time_roc_folder,'dir')
    mkdir(time_roc_folder)
end

%% Initialize stuff
ndurs = length(durations);
all_aucs = nan(2,ndurs,nb); % wake/sleep, then durations, then nbs

% Loop over wake and sleep
for iws = 1:2
    
    % Loop over durations
    for id = 1:ndurs
        fprintf('\nDoing iws %d of %d duration %d of %d\n',iws,2,id,ndurs);
        duration = durations{id};
        
        % Loop over bootstraps
        for ib = 1:nb
            
            AUC = sleep_classify_time(iws,duration);
            all_aucs(iws,id,ib) = AUC;

        end
    end
end

%% Prep duration ticks
duration_ticks = cell2mat(durations(1:ndurs-1));
new_lim = 1.2*(max(duration_ticks)-min(duration_ticks))+min(duration_ticks);
duration_ticks = [duration_ticks,new_lim];
duration_ticklabels = arrayfun(@(x) sprintf('%d',x),duration_ticks,'uniformoutput',false);
duration_ticklabels(end) = {'Full'};

if 0
figure
for iws = 1:2
    curr_aucs = squeeze(mean(all_aucs(iws,:,:),3));
    plot(duration_ticks,curr_aucs,'-','linewidth',2)
    hold on
    
end

xticks(duration_ticks)
xticklabels(duration_ticklabels)
ylabel('Mean AUC')
title('Accuracy by duration')
end
out.all_aucs = all_aucs;
out.durations = durations;
out.duration_ticklabels = duration_ticklabels;
out.duration_ticks = duration_ticks;

save([time_roc_folder,'out.mat'],'out')

end