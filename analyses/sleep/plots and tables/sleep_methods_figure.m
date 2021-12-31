function sleep_methods_figure

p = 94;% HUP197
f = 1;
time = 194084.34;
surround = 1.5;
labels = {'LA1','LA2','LA3','LD1','LD2'};
montage = 'car';
colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
name = pt(p).name;

%% Load out file and get roc stuff
%out = load([results_folder,'analysis/sleep/out.mat']);
out = load([script_folder,'analyses/sleep/data/out.mat']);
out = out.out;
disc = out.roc_out.disc;
roc = out.roc_out.roc;
auc = out.roc_out.auc;
disc_I = out.roc_out.disc_I;

%% load summ file and get ad
summ = load([results_folder,'analysis/intermediate/',name,'.mat']);
summ = summ.summ;
all_labels = summ.labels;
ekg = find_non_intracranial(all_labels);
ad = summ.ad;
ad = ad(~ekg,:);
ad = nanmean(ad,1);
ad_norm = (ad - nanmedian(ad))./iqr(ad);
wake = ad_norm > disc;
sleep = ad_norm <= disc;
all_times = summ.times;



%% Get EEG data
times = [time-surround,time+surround];
file_name = pt(p).ieeg.file(f).name;
data = download_ieeg_data(file_name,login_name,pwfile,times,1);

chLabels = data.chLabels;
values = data.values;
fs = data.fs;

%% Cleaned labels
clean_labels = decompose_labels(chLabels,name);

%% Find non-intracranial chs
non_intracranial = find_non_intracranial(clean_labels);
which_chs = find(~non_intracranial); % channels to do analysis on

%% Reject bad channels
[bad,details] = identify_bad_chs(values,which_chs,chLabels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on

%% CAR montage
[car_values,car_labels] = car_montage(values,which_chs,clean_labels);

%% Bipolar montage
[bipolar_values,~,bipolar_labels,chs_in_bipolar,which_chs_bipolar] = ...
    bipolar_montage(values,chLabels,which_chs,[],[],name);

is_run_car = ismember((1:length(clean_labels))',which_chs);
is_run_bipolar = ismember((1:length(clean_labels))',which_chs_bipolar);

%% Filters
switch montage
    case 'bipolar'
        values = bipolar_values;
        is_run = is_run_bipolar;
        curr_labels = bipolar_labels;
    case 'car'
        values = car_values;
        is_run = is_run_car;
        curr_labels = car_labels;
end

% filters
values = notch_filter(values,fs);
values = bandpass_filter(values,fs);

% make non run channels nans
run_values = values;
run_values(:,~is_run) = nan;


%% Spike detection
gdf = detector_alt(run_values,fs);

%% Match chs I want
ch_idx = ismember(clean_labels,labels);
skip = find(~ch_idx);

keep_spikes_idx = ismember(gdf(:,1),find(ch_idx));
gdf = gdf(keep_spikes_idx,:);


%% Prep figure
figure
set(gcf,'position',[10 10 1100 600])
tiledlayout(2,6,'tilespacing','tight','padding','tight')

%% Plot EEG data
nexttile([1 2])
show_spikes = 1;
show_eeg_and_spikes_select(values,gdf,fs,ch_idx,show_spikes)
xticklabels([])
yticklabels([])
title('Automated spike detection')
set(gca,'fontsize',15)

%% Plot FFT
nexttile([1 2])
ch = 40;
[ad_rat,P,freqs] = calc_ad(run_values(:,ch),fs);
freq_max = 15;
plot(freqs,P,'k-','linewidth',2)
hold on
xlim([0 freq_max])
yl = get(gca,'ylim');
ytext = yl(1) + 1.05*(yl(2)-yl(1));
new_yl = [yl(1)  yl(1) + 1.1*(yl(2)-yl(1))];
ylim = new_yl;
text(2.5,ytext,'Delta','horizontalalignment','center','fontsize',15)
text(10.5,ytext,'Alpha','horizontalalignment','center','fontsize',15)
plot([1 1],ylim,'--','color',colors(1,:),'linewidth',2)
plot([4 4],ylim,'--','color',colors(1,:),'linewidth',2)
patch([1 1 4 4],[yl(1) yl(2) yl(2) yl(1)],colors(1,:),'FaceAlpha',0.2,...
    'edgecolor',colors(1,:))

plot([8 8],ylim,'--','color',colors(2,:),'linewidth',2)
plot([13 13],ylim,'--','color',colors(2,:),'linewidth',2)
patch([8 8 13 13],[yl(1) yl(2) yl(2) yl(1)],colors(2,:),'FaceAlpha',0.2,...
    'edgecolor',colors(2,:))
xlabel('Frequency (Hz)')
ylabel('Power')
yticklabels([])
set(gca,'fontsize',15)
title('Alpha-delta power ratio calculation')

%% Alpha delta ratio over time
% Add sleep/wake
all_times = all_times/3600/24;
nexttile([1 2])
plot(all_times,ad_norm,'k-','linewidth',1)
%plot(all_times(wake),ad_norm(wake),'o','linewidth',1,'color',colors(4,:),'markersize',3)
hold on
%plot(all_times(sleep),ad_norm(sleep),'o','linewidth',1,'color',colors(3,:),'markersize',3)
plot(xlim,[disc disc],'k--','linewidth',3)
xl = get(gca,'xlim');
yl = get(gca,'ylim');
patch([xl(1) xl(2) xl(2) xl(1)],[yl(1) yl(1) disc disc],colors(4,:),'FaceAlpha',0.2)
patch([xl(1) xl(2) xl(2) xl(1)],[disc disc yl(2) yl(2)],colors(3,:),'FaceAlpha',0.2)
text(xl(2),(yl(2)+disc)/2,'Wake','rotation',270,'verticalalignment','bottom','fontsize',15,...
    'Horizontalalignment','center')
text(xl(2),(yl(1)+disc)/2,'Sleep','rotation',270,'verticalalignment','bottom','fontsize',15,...
    'Horizontalalignment','center')
xlabel('Days')
ylabel('Normalized alpha-delta ratio')
set(gca,'fontsize',15)
title('Sleep/wake classification')

%% Percent detected asleep per time of day
skip = 6;
nexttile([1 3])
all_tod_sw = out.bin_out.all_tod_sw;
tod_edges = out.bin_out.tod_edges;

nbins = size(all_tod_sw,2);
hours_mins = convert_edges_to_hours_mins(tod_edges);
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

%% ROC
nexttile([1 3])
plot(roc(:,1),roc(:,2),'k-','linewidth',2)
hold on
plot([0 1],[0 1],'k--','linewidth',2)
plot(roc(disc_I,1),roc(disc_I,2),'*','markersize',15,'linewidth',2,'color',colors(5,:));
text(roc(disc_I,1)+0.01,roc(disc_I,2)-0.05,'Sleep-wake cutoff','fontsize',15,'color',colors(5,:));
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','southeast','fontsize',15)
set(gca,'fontsize',15)
title('Classification accuracy')

print([results_folder,'analysis/sleep/methods_fig'],'-dpng')

end

function perc = convert_prop_to_perc(prop)

prop = str2num(prop);
prop = prop*100;
perc = sprintf('%d%%',prop);

end