function add_json_szs

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load json file into struct
json_file{1} = [data_folder,'sz_times/DATA_MASTER.json'];
json_file{2} = [data_folder,'sz_times/akash_sz_times/DATA_MASTER.json'];
json_file{3} = [data_folder,'sz_times/nina_sz_times/DATA_MASTER.json'];

% Loop over which json files to add
for js = 1:length(json_file)
    curr_json_file = json_file{js};

    data = loadjson(curr_json_file);
    ptnames = fieldnames(data.PATIENTS);

    for p = 1:length(ptnames)

        curr = data.PATIENTS.(ptnames{p});

        % Find the corresponding name in pt
        for j = 1:length(pt)
            if strcmp(pt(j).name,ptnames{p})

                % First, see if the sz times already available from the Erin
                % Excel source (which is preferable). If they exist, skip for
                % this reason
                nfiles = length(pt(j).ieeg.file);
                erin_done = 0;
                for f = 1:nfiles
                    if isfield(pt(j).ieeg.file(f),'sz_time_source') && ...
                            strcmp(pt(j).ieeg.file(f).sz_time_source,'Erin Excel')
                        erin_done = 1;
                    end
                end
                if erin_done
                    fprintf('\n%s has szs from Erin Excel, skipping\n',pt(j).name);
                    continue
                end

                % Second, see if there is more than one ieeg file. If there is,
                % I cannot use this, because the json file does not specify
                % which ieeg file the times are from
                
                if nfiles > 1
                    fprintf('\n%s has more than one ieeg file, skipping\n',pt(j).name);
                    continue
                end

                % Add seizure times
                if isfield(curr,'Events')
                    if isfield(curr.Events,'Ictal')
                        szs = curr.Events.Ictal;
                        szfnames = fieldnames(szs);
                        sz_times = nan(length(szfnames),2);

                        for s = 1:length(szfnames)
                            curr_sz = szs.(szfnames{s});
                            if isfield(curr_sz,'SeizureEEC') && ~isempty(curr_sz.SeizureEEC)
                                sz_start = curr_sz.SeizureEEC;
                            else
                                sz_start = curr_sz.SeizureUEO;
                            end

                            sz_end = curr_sz.SeizureEnd;
                            sz_times(s,:) = [sz_start sz_end];

                        end

                        % Skip it if sz times are funny
                        if any(sz_times == 0,'all') || any(isnan(sz_times),'all')
                            fprintf('\n%s has funny sz times, skipping\n',pt(j).name)
                            continue
                        end

                        % Fill up the sz times into the pt file
                        pt(j).ieeg.file.sz_times = sz_times;
                        pt(j).ieeg.file.sz_time_source = ['json',js];

                    end
                end

                % Quit looking for matching pt
                break
            end
        end

    end

    save([data_folder,'pt.mat'],'pt');
    
end

end