function run_to_add_coherence(whichPts)

%% Parameters
overwrite = 0;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
out_dir = [results_folder,'all_out/'];
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
            
            last_file = nan;
            last_run = nan;
            
            %% find last completed coherence
            for f = 1:length(pt(p).ieeg.file)
                exit_flag = 0;
                coherence_blocks = find(pt(p).ieeg.file(f).coherence_blocks==1);
                for ib = 1:length(coherence_blocks)
                    b = coherence_blocks(ib);
                    
                    % find first one that is not complete
                    if ~isfield(pc.file(f).run(b),'cohere_out') || isempty(pc.file(f).run(b).cohere_out)
                        last_file = f;
                        last_run = ib;
                        exit_flag = 1;
                        break
                    end
                        
                end
               
                if exit_flag == 1
                    break
                end
            end
            
            % Check if this is last one
            if isnan(last_file) && isnan(last_run)
      
                % We have finished already
                fprintf('\nAlready finished %s, skipping...\n',name)
                continue;
            end
        else
            last_file = 1;
            last_run = 1;
        end
    else
        last_file = 1;
        last_run = 1; % ok to redo one
        pc.name = name;
    end
    
    % Loop over files
    for f = last_file:length(pt(p).ieeg.file)

        file_name = pt(p).ieeg.file(f).name;
        pc.file(f).name = file_name;

        % get coherence blocks to do
        coherence_blocks = find(pt(p).ieeg.file(f).coherence_blocks==1);
        
        % Loop over cpherence blocks
        for ib = last_run:length(coherence_blocks)
            
            t = coherence_blocks(ib);
            
            tic
            fprintf('\nDoing %s file %d of %d run %d of %d\n',...
                name,f,length(pt(p).ieeg.file),ib,length(coherence_blocks));

            run_times = pt(p).ieeg.file(f).run_times(t,:);

            % Do the run
            out = individual_run_coherence(file_name,run_times,0,[],name);

            % save this
            pc.file(f).run(t).cohere_out = out;

            save([out_dir,out_name],'pc');
            ttoc = toc;
            fprintf('\nRun took %1.1f s\n',ttoc);
        end

        % Once done with all runs, set last_run to 1 before moving to next
        % file
        last_run = 1;
        
    end
    
  
end


end