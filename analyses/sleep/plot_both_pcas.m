function plot_both_pcas

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
any_pca(sleep_bins,sleep_times,locs,names)

end