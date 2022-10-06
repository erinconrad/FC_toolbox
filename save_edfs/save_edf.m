function save_edf(whichPts)

%% Parameters
overwrite = 0;
start_time = 20/24; % 8 pm first full day
duration = 12*3600; % 12 hours
chunk_size = 3600/6;
base_date = datetime('01-01-2000','InputFormat','dd-MM-yyyy');

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

    if overwrite == 1 || ~exist('meta.mat','file')
        meta.times = times;
        meta.files = cell(ntimes,1);
        meta.name = pt_name;
        meta.file_name = pt(ip).ieeg.file(1).name;
        for t = 1:ntimes
            meta.files{t} = sprintf('file%d.edf',t);
        end
        save([pt_dir,'meta.mat'],'meta');
    end

    if overwrite == 1
        last_file = 0;
    else
        listing = dir([pt_dir,'*.edf']);
        all_nums = nan(length(listing),1);
        for l = 1:length(listing)
            B = regexp(listing(l).name,'\d*','Match');
            all_nums(l) = str2num(B{1});
        end
        last_file = max(all_nums);
    end

    if last_file == ntimes
        fprintf('\nSkipping %s\n',pt_name);
    end

    if isempty(last_file)
        last_file = 0;
    end

    %% Loop over times and get data
    file_name = pt(ip).ieeg.file(1).name;
    for t = last_file+1:ntimes

        fprintf('\nDoing time %d of %d for %s\n',t,ntimes,pt_name);

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

        % Get new dates and times to update hdr
        hdr_start = seconds(times(t,1))+base_date;
        date_dashes = datestr(hdr_start,'dd-mm-yyyy');
        date_dots = datestr(hdr_start,'dd.mm.yy');
        time_dots = datestr(hdr_start,'HH.MM.SS');

        hdr.Patient = string(['filler ',date_dashes,' ',pt_name]);
        hdr.Recording = string(['Startdate ',date_dashes,' MW_1234567 MW_Inv_01 MW_Eq_01']);
        hdr.StartDate = string(date_dots);
        hdr.StartTime = string(time_dots);

        % write the edf
        edfwrite([pt_dir,meta.files{t}],hdr,values,'InputSampleType',"physical");

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

