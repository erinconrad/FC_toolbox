function ccep_pc_corr(out,cout)

%% Parameters
do_binary = 0;
plot_restricted = 0;
im = 1;
do_r2 = 1;

%% Get file locs
locations = fc_toolbox_locs;

% output folder
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'tests/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Get ccep stuff
% network and labels
ccep = cout.A;
ccep_labels_bipolar = cout.bipolar_labels;
ccep_labels_car = cout.chLabels;

% put labels into same format as pc labels
empty_labels = cellfun(@isempty,ccep_labels_bipolar);
ccep_labels_bipolar(empty_labels) = {'-'};
ccep_labels_car = cellfun(@(x) [x,'-CAR'],ccep_labels_car,'UniformOutput',false);

if im == 1
    ccep_labels = ccep_labels_bipolar;
elseif im == 2
    ccep_labels = ccep_labels_car;
end

% Identify stim and response chs
stim_chs = cout.ch_info.stim_chs;
response_chs = cout.ch_info.response_chs;

ccep_stim_labels = ccep_labels(stim_chs);
ccep_response_labels = ccep_labels(response_chs);


%% Get pc networks
pc = out.montage(im).net(1).data;
if do_r2
    pc = (pc).^2;
end
pc = wrap_or_unwrap_adjacency_fc_toolbox(pc);
pc_labels = out.montage(im).labels;
is_run = out.montage(im).is_run;
pc_labels_plot = pc_labels(is_run);


%% Restrict matrix
ccep = ccep(response_chs,stim_chs);
pc = pc(is_run,is_run);

%% Make binary
if do_binary
    
    % Binarize ccep and get percent zero
    ccep(isnan(ccep)) = 0;
    ccep(ccep~=0) = 1;
    num_zero = sum(sum(ccep==0));
    perc_zero = num_zero/(size(ccep,1)*size(ccep,2));
    
    % Make same percentage of pc zeros
    all_cors = (pc(:));
    top_prc = prctile(all_cors,perc_zero*100);
    pc((pc)<top_prc) = 0;
    pc(pc~=0) = 1;
end


%% Calculate hubness measures

% calculate indegree and outdegree
outdegree = nansum(ccep,1); % returns 1 x nchs column vector, 1 for each stim (but need to reduce chs!)
indegree = nansum(ccep,2); % returns nchs x 1 row vector
outdegree(sum(~isnan(ccep),1)==0) = nan;
indegree(sum(~isnan(ccep),2)==0) = nan;

% Make non-stim or non-response chs nans
outdegree(~stim_chs) = nan;
outdegree = outdegree';
indegree(~response_chs) = nan;

% Calculate ns
ns = nansum(pc,1);
ns(sum(~isnan(pc),1)==0) = nan;
ns = ns';
ns(~is_run) = nan;



%% Confirm labels match
if ~isequal(pc_labels,ccep_labels) 
    error('labels do not match')
end

%% Get restricted labels
restricted_x = intersect(ccep_stim_labels,pc_labels_plot);
restricted_y = intersect(ccep_response_labels,pc_labels_plot);
restricted_pc = pc(ismember(pc_labels_plot,ccep_response_labels),...
    ismember(pc_labels_plot,ccep_stim_labels));
restricted_ccep = ccep(ismember(ccep_response_labels,pc_labels_plot),...
    ismember(ccep_stim_labels,pc_labels_plot));

%% Plot
figure
set(gcf,'position',[10 10 900 900])
tiledlayout(2,2,'padding','tight','tilespacing','tight')

nexttile
if plot_restricted
    turn_nans_gray(restricted_pc)
    xticks(1:length(restricted_x))
    yticks(1:length(restricted_y))
    xticklabels(restricted_x)
    yticklabels(restricted_y)
else
    turn_nans_gray(pc)
    xticks(1:length(pc_labels_plot))
    yticks(1:length(pc_labels_plot))
    xticklabels(pc_labels_plot)
    yticklabels(pc_labels_plot)
end

nexttile
if plot_restricted
    turn_nans_gray(restricted_ccep)
    xticks(1:length(restricted_x))
    yticks(1:length(restricted_y))
    xticklabels(restricted_x)
    yticklabels(restricted_y)
else
    turn_nans_gray(ccep)
    xticks(1:length(ccep_stim_labels))
    yticks(1:length(ccep_response_labels))
    xticklabels(ccep_stim_labels)
    yticklabels(ccep_response_labels)
end

nexttile
plot(ns,outdegree,'o')
xlabel('NS')
ylabel('CCEPs outdegree')
[r,p] = corr(ns,outdegree,'rows','pairwise');
pause(0.2)
yl = ylim;
xl = xlim;
pause(0.2)
text(xl(1),yl(2),sprintf('r = %1.2f, p = %1.3f',r,p),'verticalalignment','top')

nexttile
plot(ns,indegree,'o')
xlabel('NS')
ylabel('CCEPs indegree')
[r,p] = corr(ns,indegree,'rows','pairwise');
pause(0.2)
yl = ylim;
xl = xlim;
pause(0.2)
text(xl(1),yl(2),sprintf('r = %1.2f, p = %1.3f',r,p),'verticalalignment','top')




end