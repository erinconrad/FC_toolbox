function fc_ccep_corr_plot(all)

%% Parameters
do_bin = 0;
do_log = 0;
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
pc = all.net.pc.data;
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

%% Get outdegree, indegree, ns
outdegree = (nansum(ccep,1))';
indegree = nansum(ccep,2);
ns = nansum(pc,2);

% sanity check
if abs(nansum(pc,2)-(nansum(pc,1))')>1e-6, error('oh no'); end

if strcmp(corr_type,'Spearman')
    rtext = '\rho';
else
    rtext = 'r';
end

%% Scatter plots
figure
set(gcf,'position',[10 10 1000 800])
tiledlayout(2,2,'tilespacing','compact','padding','tight')

nexttile
turn_nans_gray(ccep)
xticklabels([])
yticklabels([])
xlabel('Stim')
ylabel('Response')
if do_bin
    title('CCEP (binarized)')
elseif do_log
    title('CCEP (log)')
else
    title('CCEP');
end
set(gca,'fontsize',15)

nexttile
turn_nans_gray(pc)
xticklabels([])
yticklabels([])
title('Time-averaged pearson correlation')
set(gca,'fontsize',15)

nexttile

if show_labels
    plot(ns,indegree,'o','linewidth',2,'markersize',10,'color','w')
    hold on
    text(ns,indegree,labels,'horizontalalignment','center',...
        'verticalalignment','middle','fontsize',15)
else
    plot(ns,indegree,'o','linewidth',2,'markersize',10)
end
xlabel('PC node strength')
ylabel('CCEP indegree');
set(gca,'fontsize',15)
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('%s = %1.2f',rtext,corr(ns,indegree,'Type',corr_type)),...
    'verticalalignment','top','fontsize',20,'color','k')

nexttile
if show_labels
    plot(ns,outdegree,'o','linewidth',2,'markersize',10,'color','w')
    hold on
    text(ns,outdegree,labels,'horizontalalignment','center',...
        'verticalalignment','middle','fontsize',15)
else
    plot(ns,outdegree,'o','linewidth',2,'markersize',10)
end
xlabel('PC node strength')
ylabel('CCEP outdegree');
set(gca,'fontsize',15)
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('%s = %1.2f',rtext,corr(ns,outdegree,'Type',corr_type)),...
    'verticalalignment','top','fontsize',20,'color','k')


print([out_folder,pt_name,bin_text],'-dpng');
end