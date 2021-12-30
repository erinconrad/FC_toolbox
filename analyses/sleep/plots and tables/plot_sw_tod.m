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
set(gcf,'position',[440 296 905 500])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

%% Proportion asleep by time of day
nexttile
all_tod_sw = bin_out.all_tod_sw;
tod_edges = bin_out.tod_edges;

nbins = size(all_tod_sw,2);
hours_mins = convert_edges_to_hours_mins(tod_edges);
tod_sw = squeeze(nansum(all_tod_sw,1));
ind_pt_prop = all_tod_sw(:,:,2)./(all_tod_sw(:,:,2)+all_tod_sw(:,:,1));
prop_asleep = squeeze(nanmean(ind_pt_prop,1));

polar = convert_times_to_polar(tod_edges,'radians');
polarhistogram('BinEdges',polar,'BinCounts',prop_asleep,...
    'edgecolor','none','facealpha',0.5)%'displaystyle','stairs')
set(gca,'ThetaDir','clockwise');
set(gca,'ThetaZeroLocation','top');
rlabs = get(gca,'rticklabel');
rlabs = cellfun(@convert_prop_to_perc,rlabs,'uniformoutput',false);
set(gca,'rticklabel',rlabs)
thetaticks(polar(1:skip:nbins)*360/(2*pi))
thetaticklabels(hours_mins(1:skip:nbins+1))
set(gca,'fontsize',15)
title('Percent detected asleep')

%{
plot(prop_asleep,'k','linewidth',2)
xticks(1:skip:nbins)
xticklabels(hours_mins(1:skip:nbins+1))
ylabel('Proportion detected asleep')
xlabel('Time of day')
set(gca,'fontsize',15)
title('Proportion detected asleep by time of day')
%}


%% prep stuff
loc = circ_out.all_locs;
all_tod_rate = circ_out.all_tod_rate;
norm_rate = (all_tod_rate - nanmedian(all_tod_rate,2))./iqr(all_tod_rate,2);
iqr_tod_rate = [prctile(norm_rate,25,1);prctile(norm_rate,75,1)];
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');


median_temporal = nanmedian(norm_rate(temporal,:),1);
iqr_temporal = [prctile(norm_rate(temporal,:),25,1);prctile(norm_rate(temporal,:),75,1)];
median_extra = nanmedian(norm_rate(extra,:),1);
iqr_extra = [prctile(norm_rate(extra,:),25,1);prctile(norm_rate(extra,:),75,1)];

%% Spike rate by time of day
nexttile
median_tod_rate = (nanmedian(norm_rate,1));

polarhistogram('BinEdges',polar,'BinCounts',median_tod_rate+min(median_tod_rate)+1,...
    'edgecolor','none')
set(gca,'ThetaDir','clockwise');
set(gca,'ThetaZeroLocation','top');
thetaticks(polar(1:skip:nbins)*360/(2*pi))
thetaticklabels(hours_mins(1:skip:nbins+1))
set(gca,'fontsize',15)
title('Spike rate')


%% Spike rate temporal
nexttile([1 2])
%shaded_error_bars(1:length(median_temporal),median_temporal,iqr_temporal,myColours(1,:));
plot(median_temporal,'linewidth',2);
hold on
plot(median_extra,'linewidth',2);
xticks(1:skip:nbins)
xticklabels(hours_mins(1:skip:nbins+1))
ylabel('Normalized spike rate')
xlabel('Time of day')
legend({'Temporal','Extra-temporal'},'fontsize',15)
set(gca,'fontsize',15)
title('Spike rate by seizure localization')




%{
shaded_error_bars(1:length(median_tod_rate),median_tod_rate,iqr_tod_rate,'k');
plot(tod_rate,'linewidth',2)
xticks(1:skip:nbins)
xticklabels(hours_mins(1:skip:end))
ylabel('Normalized spike rate')
xlabel('Time of day')
set(gca,'fontsize',15)
title('Spike rate by time of day')
%}




%{
counts = nansum(circ_out.all_tod_rate,1);
observations = convert_counts_to_observations(counts,tod_edges);
polar = convert_times_to_polar(observations,'radians');
circ_plot(polar,'hist',[],length(tod_edges),true,true,'linewidth',2,'color','r')
[pval z] = circ_rtest(polar);

% Make fake random data
a = 0;
b = 2*pi;
fake_observations = (b-a).*rand(length(polar),1)+a;
[pval z] = circ_rtest(fake_observations);
%}
end

function perc = convert_prop_to_perc(prop)

prop = str2num(prop);
prop = prop*100;
perc = sprintf('%d%%',prop);

end