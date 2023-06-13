function mt_spike_validation(whichPts,overwrite)


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_summ_path = [results_folder,'edf_summ_out/'];
edf_path = [results_folder,'edf_out/'];
data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load pt folder
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% LOad validation file
validation_file = [scripts_folder,'spike_detector/Manual validation.xlsx'];
mT = readtable(validation_file,'Sheet','strange_elec_names');

if isempty(whichPts)
    whichPts = 1:length(pt);
end

for i = 1:length(whichPts)
    ip = whichPts(i);
    name = pt(ip).name;

    % see if I've already done it
    if exist([edf_summ_path,name,'/bipolar.png'],'file') ~=0
        if overwrite == 1
            fprintf('\nOverwriting %s\n',name);
        else
            fprintf('\nSkipping %s\n',name);
            continue
        end
    end

    % Load the summary file
    if exist([edf_summ_path,name,'/summ.mat'],'file') == 0
        fprintf('\nNo summary file for %s, skipping\n',name);
        continue
    end
    out = load([edf_summ_path,name,'/summ.mat']);
    out = out.out;

    all_spike_times = out.all_spike_times;
    labels = out.labels;
    montages = out.montages;
    nmontages = length(montages);

    % Plot random spike detections
    for im = 1:nmontages
        %plot_random_spikes(all_spike_times{im},name,labels,montages{im},edf_path,edf_summ_path)
            plot_random_spikes(all_spike_times{im},name,out.labels,montages{im},...
                edf_path,edf_summ_path,mT,pt,ip,0)

    end
end


end