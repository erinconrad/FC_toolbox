function plot_sw_tod

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


skip = 6;

figure

nexttile
all_tod_sw = bin_out.all_tod_sw;
tod_edges = bin_out.tod_edges;
hours_mins = convert_edges_to_hours_mins(tod_edges);
tod_sw = squeeze(nansum(all_tod_sw,1));
ind_pt_prop = all_tod_sw(:,:,2)./(all_tod_sw(:,:,2)+all_tod_sw(:,:,1));
%prop_asleep = tod_sw(:,2)./(tod_sw(:,2) + tod_sw(:,1));
prop_asleep = squeeze(nanmean(ind_pt_prop,1));

plot(prop_asleep,'linewidth',2)
xticks(1:skip:size(tod_sw,1))
xticklabels(hours_mins(1:skip:size(tod_sw,1)))
ylabel('Proportion asleep')
xlabel('Time of day')
set(gca,'fontsize',15)
title('Proportion asleep by time of day')

nexttile
all_tod_rate = circ_out.all_tod_rate;
norm_rate = (all_tod_rate - nanmedian(all_tod_rate,2))./iqr(all_tod_rate,2);
tod_rate = squeeze(nanmedian(norm_rate));
plot(tod_rate,'linewidth',2)
xticks(1:skip:length(tod_rate))
xticklabels(hours_mins(1:skip:length(tod_rate)))
ylabel('Spike rate')
xlabel('Time of day')
set(gca,'fontsize',15)
title('Spike rate by time of day')


nexttile
eleven_to_five = circ_out.eleven_to_five;
plot_paired_data(eleven_to_five',{'Night','Day'},'Spike rate','paired',plot_type)
%{
nexttile
nine_pm_to_nine_am = tod_edges(2:end) < 6*3600 | tod_edges(2:end) > 23*3600;
early_late_rate = [nanmean(all_tod_rate(:,nine_pm_to_nine_am),2),...
    nanmean(all_tod_rate(:,~nine_pm_to_nine_am),2)];
plot_paired_data(early_late_rate',{'Night','Day'},'Spike rate','paired',plot_type)
%}

figure
%{
asleep = ind_pt_prop;
asleep = (asleep - nanmean(asleep,2))./nanstd(asleep,[],2);
[coeff,score,latent] = pca(asleep);
%}

rate = all_tod_rate;
rate = (rate - nanmean(rate,2))./nanstd(rate,[],2);
[coeff,score,latent] = pca(rate);

nexttile
plot(coeff(:,1),'linewidth',2)
xticks(1:skip:length(tod_rate))
xticklabels(hours_mins(1:skip:length(tod_rate)))
ylabel('Component 1')

nexttile
plot(coeff(:,2),'linewidth',2)
xticks(1:skip:length(tod_rate))
xticklabels(hours_mins(1:skip:length(tod_rate)))
ylabel('Component 2')






end