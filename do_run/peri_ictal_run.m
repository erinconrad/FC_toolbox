function peri_ictal_run(whichPts)

%% Parameters
overwrite = 0;
tw = 2; % 2 second time window for pc calculations
which_net = 'pc';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
out_dir = [results_folder,'peri_ictal_out/'];
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

if isempty(whichPts)
    whichPts = 1:length(pt);
end

if ischar(whichPts)
    for j = 1:length(pt)
        if strcmp(whichPts,pt(j).name)
            whichPts = j;
            break
        end
    end
end

if iscell(whichPts)
    whichPtsNames = whichPts;
    whichPts = nan(length(whichPts));
    for i = 1:length(whichPts)
        for j = 1:length(pt)
            if strcmp(whichPtsNames{i},pt(j).name)
                whichPts(i) = j;
                break
            end
        end
    end
end


%% Do the run
% Loop over pts
for i = 1:length(whichPts)
    
    clear pc % so it doesn't have old stuff
    
    p = whichPts(i);
    name = pt(p).name;
    out_name = [name,'_pc.mat'];
    
    % See if it already exists
    if overwrite == 0
        if exist([out_dir,out_name],'file') ~= 0
            pc = load([out_dir,out_name]);
            pc = pc.pc;
            
            last_file = length(pc.file);
            last_sz = length(pc.file(last_file).sz);
            last_run = length(pc.file(last_file).sz(last_sz).run);
            
            if last_file == length(pt(p).ieeg.file) && ...
                    last_sz == length(pt(p).ieeg.file(last_file).sz_run_times) && ...
                    last_run == size(pt(p).ieeg.file(last_file).sz_run_times{last_sz},1)
                % We have finished already
                fprintf('\nAlready finished %s, skipping...\n',name)
                continue;
            end
        else
            last_file = 1;
            last_sz = 1;
            last_run = 1;
            pc.name = name;
        end
    else
        last_file = 1;
        last_sz = 1;
        last_run = 1; % ok to redo one
        pc.name = name;
    end
    
    % Loop over files
    for f = last_file:length(pt(p).ieeg.file)

        file_name = pt(p).ieeg.file(f).name;
        pc.file(f).name = file_name;

        % Loop over seizures within file
        for s = last_sz:length(pt(p).ieeg.file(f).sz_run_times)
            
            % loop over run times within seizure
            for t = last_run:size(pt(p).ieeg.file(f).sz_run_times{s},1)
            
                tic
                fprintf('\nDoing %s file %d of %d seizure %d of %d run %d of %d\n',...
                    name,f,length(pt(p).ieeg.file),...
                    s,length(pt(p).ieeg.file(f).sz_run_times),...
                    t,size(pt(p).ieeg.file(f).sz_run_times{s},1));

                run_times = pt(p).ieeg.file(f).sz_run_times{s}(t,:);
                sz_times = pt(p).ieeg.file(f).sz_times(s,:);

                % Do the run
                out = individual_run(file_name,run_times,tw,which_net,0,[],name);

                % save this
                pc.file(f).sz(s).run(t).data = out;
                pc.file(f).sz(s).run(t).run_times = run_times;
                pc.file(f).sz(s).sz_times = sz_times;

                save([out_dir,out_name],'pc');
                ttoc = toc;
                fprintf('\nRun took %1.1f s\n',ttoc);
            
            end
            
            % Once done with all runs, set last_run to 1 before moving to next
            % file
            last_run = 1;
        end

        % Once done with all szs, set last_sz to 1 before moving to next
        % file
        last_sz = 1;
        
    end
    % Once done with all files, set last_file to 1 before moving to next
    % patient
    last_file = 1;
    
end


end