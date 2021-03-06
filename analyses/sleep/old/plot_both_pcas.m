function plot_both_pcas

%% Parameters

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

%% Load the out file
load([out_folder,'out.mat']);

%% Get stuff for both pcas
sleep_bins = out.sleep_hist_out.all_pts_spikes_bins;
sleep_times = out.sleep_hist_out.times;
locs = out.circ_out.all_locs;
names = out.circ_out.names;
sz_bins = out.sz_out.all_pts_spikes_bins;
sz_times = out.sz_out.times;

%% Do both pcas
any_pca(sleep_bins,sleep_times,locs,names,'Spike patterns around sleep')
print([out_folder,'sleep_pca'],'-dpng')
close(gcf)
any_pca(sz_bins,sz_times,locs,names,'Peri-ictal spike patterns')
print([out_folder,'seizure_pca'],'-dpng')
close(gcf)

%% Also show all histograms
%plot_all_histograms(sleep_bins,sleep_times,names)
%plot_all_histograms(sz_bins,sz_times,names)

end