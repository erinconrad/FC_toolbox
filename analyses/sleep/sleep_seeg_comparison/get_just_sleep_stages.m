function get_just_sleep_stages

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_folder = [results_folder,'edf_out/'];
out_folder = [results_folder,'sleepseeg_stages/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

listing = dir(edf_folder);
for l = 1:length(listing)
    if exist([edf_folder,listing(l).name,'/sleep_stage.mat'])
        sout = load([edf_folder,listing(l).name,'/sleep_stage.mat']);
        sout = sout.sout;

        % Get times of sleep transitions
        st_dates= sout.Summary(2:end,2);
        st_times = sout.Summary(2:end,3);
        stage = sout.Summary(2:end,4);
        seeg_secs = convert_transitions_to_ieeg_secs(st_dates,st_times);

        out.name = listing(l).name;
        out.ChL = sout.ChL;
        out.Summary = sout.Summary;
        out.SleepStage = sout.SleepStage;
        out.stage = stage;
        out.st_dates = st_dates;
        out.st_times = st_times;
        out.seeg_secs = seeg_secs;
        out.which_file = 1;
        out.duration_hours = 12;

        save([out_folder,out.name,'.mat'],'out');
    end
end

end
