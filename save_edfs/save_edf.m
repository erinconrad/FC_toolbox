function save_edf(whichPts)

%% Parameters
start_time = 20/24; % 8 pm first full day
duration = 12*3600; % 12 hours
chunk_size = 3600/6;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
out_dir = [results_folder,'edf_out/'];
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

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

for ip = whichPts
    
    %% Decide start time

    % Get start time of day of first file
    ff_start = pt(ip).ieeg.file(1).start_time;

    % download start should be 8 pm the first full day
    dl_start = (start_time-ff_start)*3600*24 + 3600*24;

    %% Prepare and save metadata
    pt_name = pt(ip).name;
    times = nan(duration/chunk_size,2);
    ntimes = size(times,1);
    for t = 1:ntimes
        times(t,1) = dl_start + (t-1)*chunk_size;
        times(t,2) = dl_start+t*chunk_size;
    end
    
    pt_dir = [out_dir,pt_name,'/'];
    if ~exist(pt_dir,'dir')
        mkdir(pt_dir)
    end
    meta.times = times;
    meta.files = cell(ntimes,1);
    meta.name = pt_name;
    meta.file_name = pt(ip).ieeg.file(1).name;
    for t = 1:ntimes
        meta.files{t} = sprintf('file%d.edf',t);
    end
    save([pt_dir,'meta.mat'],'meta');

    %% Loop over times and get data
    file_name = pt(ip).ieeg.file(1).name;
    for t = 1:ntimes
        % Get ieeg data
        data = download_ieeg_data(file_name,login_name,pwfile,times(t,:),1);
        values = data.values;
        values = floor(values);
        fs = data.fs;
        

        % Downsample the data
        desired_fs = 256;
        n = fs/desired_fs;
        values = downsample(values,n);
        
        labels = decompose_labels(data.chLabels(:,1),pt_name);
        ekg = find_non_intracranial(labels);

        nchs = sum(~ekg);
        nsamples = size(values,1);
        values = values(:,~ekg);
        labels = labels(~ekg);
       
        %values = values/1000; % go from uV to mV
         
        % Prep header
        hdr = edfheader("EDF+");
        hdr.NumDataRecords = 1;
        hdr.DataRecordDuration = seconds(nsamples/desired_fs);
        hdr.NumSignals = nchs;
        hdr.SignalLabels = labels;
        hdr.PhysicalDimensions = repelem("uV",nchs);
        hdr.PhysicalMin = min(values);
        hdr.PhysicalMax = max(values);
        hdr.DigitalMin = repmat(-32768,1,nchs);
        hdr.DigitalMax = repmat(32767,1,nchs);

        % write the edf
        edfw = edfwrite([pt_dir,meta.files{t}],hdr,values,'InputSampleType',"physical");

        % Update meta file
        if t == 1
            meta.ekg = ekg;
            meta.labels = labels;
            meta.original_fs = fs;
            meta.new_fs = desired_fs;
            save([pt_dir,'meta.mat'],'meta');
        end

    end

end

end

