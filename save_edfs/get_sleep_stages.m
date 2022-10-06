function get_sleep_stages

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
out_dir = [results_folder,'edf_out/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Add sleep seeg folder to path
addpath(genpath(locations.sleep_seeg_folder))

% Loop over folders in edf out
listing = dir(out_dir);
for l = 1:length(listing)
    
    % loop over edf files within patient-specific folder
    all_edf_files = {};
    sub_listing = dir([out_dir,listing(l).name,'/*.edf']);
    for f = 1:length(sub_listing)
        all_edf_files = [all_edf_files;[out_dir,listing(l).name,'/',sub_listing(f).name]];
    end

    if isempty(all_edf_files)
        continue
    end

    all_nums = nan(length(all_edf_files),1);
    for i = 1:length(all_edf_files)
        B = regexp(all_edf_files{i},'\d*','Match');
        all_nums(i) = str2num(B{2});
    end
    [~,I] = sort(all_nums);
    all_edf_files = all_edf_files(I);

    % Run the sleep SEEG
    [Summary,SleepStage,ChL]=SleepSEEG_erin(all_edf_files(1:10),0);

    sout.Summary = Summary;
    sout.SleepStage = SleepStage;
    sout.ChL = ChL;

    % save the output
    save([out_dir,listing(l).name,'/sleep_stage.mat'],'sout')

end


end