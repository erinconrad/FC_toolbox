function edf_visualize_data(name,file)

% I plotted a few random times and it seems to align with expected

%% params
which_ch = 4;
secs_to_plot = 15;
ref_start_time = datetime('01/01/2000 00:00:00','InputFormat','MM/dd/yyyy hh:mm:ss');

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
out_dir = [results_folder,'edf_out/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Load the edf file
file_path = [out_dir,name,'/file',sprintf('%d',file),'.edf'];
info = edfinfo(file_path);
tt = edfread(file_path);


%% Get time on ieeg.org in seconds
start_time = datetime([char(info.StartDate),' ',char(info.StartTime)],...
    'InputFormat','dd.MM.yy HH.mm.ss');

td = start_time-ref_start_time;
tds = seconds(td);

%% Plot the signal on the first channel
fs = info.NumSamples/seconds(info.DataRecordDuration);
recnum = 1;
signum = which_ch;
samples_to_plot = round(secs_to_plot/seconds(info.DataRecordDuration)*info.NumSamples(signum));

t = (0:samples_to_plot)/fs(signum);
y = tt.(signum){recnum}(1:samples_to_plot+1);

plot(t,y)
legend(strcat("Record ",int2str(recnum),", Signal ",info.SignalLabels(signum)))
title(sprintf('%d seconds',tds))
file_path
fprintf('\n%d seconds\n',tds)

end