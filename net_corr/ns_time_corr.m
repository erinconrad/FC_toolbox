function ns_time_corr(all)

%% Parameters
do_bin = 0;
do_log = 1;
corr_type = 'Pearson';
show_labels = 0;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'corrs/fc_ccep/'];
script_folder = locations.script_folder;

addpath(genpath(script_folder))

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

pt_name = all.pt_name;
labels = all.labels;

%% Get two networks
ns = all.ns_time;
nruns = all.nruns;
block_dur = all.block_dur_secs;
end_time = block_dur*nruns/3600/24;
ccep = all.net.ccep.data;

%% Binarize ccep
if do_bin
    a = ~isnan(ccep);
    ccep(a) = 1;
    ccep(~a) = 0;
    
    bin_text = '_bin';
elseif do_log
    ccep = log(ccep);
    bin_text = '_log';
else
    bin_text = '';
end

%% Get outdegree, indegree
outdegree = (nansum(ccep,1))';
indegree = nansum(ccep,2);

%% Correlation over time
out_ns_time = corr(outdegree,ns,'Type',corr_type,'rows','pairwise');
in_ns_time = corr(indegree,ns,'Type',corr_type,'rows','pairwise');

if strcmp(corr_type,'Spearman')
    rtext = '\rho';
else
    rtext = 'r';
end

figure
tiledlayout(2,1,'tilespacing','compact','padding','tight')

nexttile()
h = turn_nans_gray(ns);
set(h,'XData',[0:end_time]);
%xticks(1:floor(end_time))
xticklabels([])
xlim([0 end_time])
%xlabel('Days')
eo = 1:2:length(labels);
yticks(eo)
yticklabels(labels(eo))
set(gca,'fontsize',15)
title('Pearson correlation node strength')
colorbar

%{
nexttile
plot(out_ns_time)
xlabel('Time')
ylabel('Outdegree-NS correlation')
%}

nexttile
plot(linspace(0,end_time,length(in_ns_time)),in_ns_time)
xlim([0 end_time])
xticks(1:floor(end_time))
xlabel('Days')
ylabel('Indegree-NS correlation')
set(gca,'fontsize',15)

print(gcf,[out_folder,pt_name,'_time'],'-dpng')

end