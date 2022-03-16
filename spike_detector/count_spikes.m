function count_spikes

%% Locations
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
spike_folder = [results_folder,'all_out/'];
data_folder = [locations.main_folder,'data/'];
out_folder = [results_folder,'sp_validation/'];
out_file = [out_folder,'spike_summary.xlsx'];

%% get table of names
nameT = readtable(out_file);

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;


whichPts = [];
listing = dir([spike_folder,'*.mat']);
for i = 1:length(listing)
    C = listing(i).name;
    temp_name = strsplit(C,'_');
    temp_name = temp_name{1};
    for j = 1:length(pt)
        pt_name = pt(j).name;
        if strcmp(temp_name,pt_name)
            whichPts = [whichPts,j];
            break
        end
    end
end

for p = whichPts
    pt_name = pt(p).name;
    
    %% Load spike file
    out = load([spike_folder,sprintf('%s_pc.mat',pt_name)]);
    out = out.pc;
    
    %% concatenate all spikes into one long thing
    % Include an extra column for the file index and block
    for im = 1:2
        all_spikes = [];
        
        if ~isfield(out.file(1).run(1),'run_times')
            nblocks = nan;
            nspikes = nan;
            spikes_per_min = nan;
        else
        
            mins_per_block = diff(out.file(1).run(1).run_times)/60; % divide by 60 to convert from secs to mins
            nblocks = 0;
            ndup = 0;
            for f = 1:length(out.file)

                for h = 1:length(out.file(f).run)
                    gdf = out.file(f).run(h).data.montage(im).spikes;

                    % remove duplicates
                    if ~isempty(gdf)
                        [gdf,n_removed] = remove_duplicates(gdf);
                        ndup = ndup + n_removed;
                    end

                    nblocks = nblocks + 1;
                    all_spikes = [all_spikes;gdf,...
                        repmat(f,size(gdf,1),1),...
                        repmat(h,size(gdf,1),1)];
                end
            end

            if ndup >0, error('why'); end

            nspikes = size(all_spikes,1);
            spikes_per_min = nspikes/nblocks/mins_per_block;
        end

        %% Prep table
        if im == 1
            T = table(nblocks,nspikes,spikes_per_min);
        elseif im == 2
            T = table(nspikes,spikes_per_min);
        end
        %% write to csv file
        row = p+1;

        % confirm the row is as expected
        if ~strcmp(nameT.name{p},pt_name), error('why'); end

        % Get names of rows and columns to print
        if im == 1
            print_range = sprintf('B%d:D%d',row,row);
        elseif im == 2
            print_range = sprintf('G%d:H%d',row,row);
        end

        % write the info
        writetable(T,out_file,'Range',print_range,'WriteVariableNames',false);
    end
    
end


end