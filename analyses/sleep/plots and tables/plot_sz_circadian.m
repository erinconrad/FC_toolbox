function plot_sz_circadian

%% Parameters
plot_type = 'scatter';
nblocks = 6;
myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];



locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

%% Load out file and get roc stuff
out = load([out_folder,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];


figure

nexttile
periods = sz_circ_out.periods;
iqr_psd = sz_circ_out.iqr_psd;
all_psd = sz_circ_out.all_psd;
median_psd = nanmedian(all_psd,1);
mp = shaded_error_bars(periods,median_psd,iqr_psd,[0 0 0]);

nexttile
pre_wake = sz_circ_out.pre_wake;
for i = 1:length(pre_wake)
    plot(i,sum(pre_wake{i}==1)/(sum(pre_wake{i}==1) + sum(pre_wake{i}==0)),'o')
    hold on
end

end