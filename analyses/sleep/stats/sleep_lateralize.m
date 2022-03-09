function sleep_lateralize

%% Parameters
min_num_elecs_each_side = 5;

%% File lcos
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
%out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

rate_sw = bin_out.all_elecs_rates_sw;
labels = bin_out.all_elecs_names;
epi_lat = circ_out.all_lats;

%% Get pairs of L-R spike rates in sleep
npts = length(rate_sw);
rate_sleep_lr = nan(npts,2);
for ip = 1:npts
    
    % make sure all intracranial
    curr_labels = labels{ip};
    assert(sum(find_non_intracranial(curr_labels)) == 0)
    
    % Parse the laterality
    curr_lats = parse_electrode_laterality(curr_labels);
    
    % Skip if fewer than 5 for either laterality
    if sum(strcmp(curr_lats,'L')) < min_num_elecs_each_side || ...
            sum(strcmp(curr_lats,'R')) < min_num_elecs_each_side
        continue;
    end
    
    % Get current spike rates in sleep of all elecs
    curr_sleep_rates = rate_sw{ip}(:,2);
    
    % Spike rate in sleep for left sided and right side
    rate_sleep_lr(ip,1) = nanmean(curr_sleep_rates(strcmp(curr_lats,'L')));
    rate_sleep_lr(ip,2) = nanmean(curr_sleep_rates(strcmp(curr_lats,'R')));
    
end

%% Remove those with nans
any_nan = any(isnan(rate_sleep_lr),2);
rate_sleep_lr(any_nan,:) = [];
epi_lat(any_nan) = [];

%% Remove those that are bilateral
is_bilat = ~strcmp(epi_lat,'left') & ~strcmp(epi_lat,'right');
epi_lat(is_bilat) = [];
rate_sleep_lr(is_bilat,:) = [];

if 0
    table(epi_lat,rate_sleep_lr)
end

%% What % of patients has higher spike rate on side of soz?
higher_on_left = rate_sleep_lr(:,1) > rate_sleep_lr(:,2);
same_side = (higher_on_left & strcmp(epi_lat,'left')) | ...
    (~higher_on_left & strcmp(epi_lat,'right'));

fprintf('\n%d of %d (%1.1f%%) patients had more spikes in sleep on the side of the SOZ.\n',...
    sum(same_side),length(same_side),sum(same_side)/length(same_side)*100);

end