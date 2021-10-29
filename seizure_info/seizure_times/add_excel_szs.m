function add_excel_szs

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Get sheet names from sz table
sheet = sheetnames([data_folder,'sz_times/more_seizure_times.xlsx']);


% Loop over patients
for p = 1:length(pt)
    name = pt(p).name;
    
    for s = 1:length(sheet)
        
        sname = sheet(s);
        if contains(sname,name)
            if ~contains(sname,'complete')
                break
            else
                
                % Load sheet
                T = readtable([data_folder,'sz_times/more_seizure_times.xlsx'],'Sheet',sname);
                
                % initialize data
                for f = 1:length(pt(p).ieeg.file)
                    pt(p).ieeg.file(f).sz_info = [];
                end
                
                % fill sz info
                for sz = 1:size(T,1)
                    f = T.IEEG_file(sz);
                    info.eec = T.EEC(sz);
                    info.ueo = T.UEO(sz);
                    info.end = T.End(sz);
                    info.elec = T.OnsetElectrode(sz);
                    info.type = T.Type(sz);
                    info.notes = T.Notes(sz);
                    
                    pt(p).ieeg.file(f).sz_info = [pt(p).ieeg.file(f).sz_info;
                        info];
                end
                
                % Get all times in file
                for f = 1:length(pt(p).ieeg.file)
                    sz_times = nan(length(pt(p).ieeg.file(f).sz_info),2);
                    for z = 1:length(pt(p).ieeg.file(f).sz_info)
                        eec = pt(p).ieeg.file(f).sz_info(z).eec;
                        ueo = pt(p).ieeg.file(f).sz_info(z).ueo;
                        send = pt(p).ieeg.file(f).sz_info(z).end;
                        
                        % preferentially use eec
                        if isnan(eec)
                            sz_times(z,:) = [ueo send];
                        else
                            sz_times(z,:) = [eec send];
                        end
                    end
                    pt(p).ieeg.file(f).sz_times = sz_times;
                    
                    % Note the source of the seizures
                    pt(p).ieeg.file(f).sz_time_source = 'Erin Excel';
                end
                
                % Quit looking for matching pt
                break
            
            end
            
            
        end


    end
    
end

save([data_folder,'pt.mat'],'pt');


end