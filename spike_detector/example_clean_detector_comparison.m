%% Example clean_detector_fast

clear; clf; clc;

%%%% Load in an example data set
%  From
%       https://rdcu.be/eukIo
%       Lin et al. An EEG dataset for interictal epileptiform discharge with spatial distribution information.
%       Sci Data 12, 229 (2025). https://doi.org/10.1038/s41597-025-04572-1
%  Download data
%       https://doi.org/10.6084/m9.figshare.28069568

load('/Users/mak/Downloads/opensource-dataset/MAT_Files/DA00100B.mat')
Fs  = 500;

%%%% Scale to uV for detector.
eeg_data = eeg_data*1e6;

%%%% Choose one channel to analyze.
channel_to_analyze = 1;

%%%% Run new detector
d0  = eeg_data(channel_to_analyze,:);
tic
gdf = clean_detector_one_channel_fast(transpose(d0),Fs);
elapsedTimeNew = toc;

%%%% Run original detector.
tic
gdf_original = clean_detector(transpose(eeg_data),Fs);
elapsedTimeOrig  = toc;
elapsedTimeOrig_per_channel = elapsedTimeOrig/size(eeg_data,1);

%%%% Plot diagnostics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% New detector.
dspikes = [];
for k=1:size(gdf,1)
    start_time = gdf(k,2)-Fs/2;
    stop_time  = gdf(k,2)+Fs/2-1;
    dspikes(k,:) = d0(start_time:stop_time);
end
taxis = (-Fs/2:Fs/2-1)/Fs;
subplot(2,2,1)
plot(taxis, dspikes')
title(['New Method. Counts ' num2str(size(gdf,1)) ', Time ' num2str(elapsedTimeNew,2) ' s'])
grid on

subplot(2,2,3)
med   = median(dspikes,1);
[r,q] = iqr(dspikes,1);
plot(taxis, med, 'k', 'LineWidth', 2);
hold on
plot(taxis, q, '--r')
hold off
grid on

%%%% Original detector.
indices = find(gdf_original(:,1)==channel_to_analyze);
gdf_original = gdf_original(indices,:);
dspikes = [];
for k=1:size(gdf_original,1)
    start_time = gdf_original(k,2)-Fs/2;
    stop_time  = gdf_original(k,2)+Fs/2-1;
    dspikes(k,:) = d0(start_time:stop_time);
end
taxis = (-Fs/2:Fs/2-1)/Fs;
subplot(2,2,2)
plot(taxis, dspikes')
title(['Original Method. Counts ' num2str(size(gdf_original,1)) ', Time ' num2str(elapsedTimeOrig_per_channel,2) ' s per electrode'])
grid on

subplot(2,2,4)
med   = median(dspikes,1);
[r,q] = iqr(dspikes,1);
plot(taxis, med, 'k', 'LineWidth', 2);
hold on
plot(taxis, q, '--r')
hold off
grid on
