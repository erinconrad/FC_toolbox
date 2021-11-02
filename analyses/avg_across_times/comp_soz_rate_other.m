function comp_soz_rate_other

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/avg_spikes_over_time/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the out file
out = load([out_folder,'out.mat']);
out = out.out;

npts = length(out);
missing_soz =  zeros(npts,1);
spikes_soz_not = nan(npts,2);

for p = 1:length(out)
    
    soz_labels = out(p).soz_labels;
    labels = out(p).labels;
    avg_spikes = out(p).avg_spikes;
    
    if isempty(soz_labels)
        missing_soz(p) = 1;
        continue
    end
    
    % get indices of soz labels
    is_soz = ismember(labels,soz_labels);
    spikes_soz = nanmean(avg_spikes(is_soz));
    spikes_not = nanmean(avg_spikes(~is_soz));
    
    spikes_soz_not(p,:) = [spikes_soz spikes_not];
    
end

% Plot
 plot_paired_data(spikes_soz_not',{'SOZ','Not'},...
     'Spikes/elec/min','paired','errorbar')
 print(gcf,[out_folder,'soz_comp'],'-dpng')
 

end