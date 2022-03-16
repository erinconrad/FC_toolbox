function all_pt_spike_rates

%% Locations
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
data_folder = [locations.main_folder,'data/'];

spike_folder = [results_folder,'all_out/'];
out_folder = [results_folder,'sp_validation/'];


if ~exist(out_folder,'dir'), mkdir(out_folder); end

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

%% Prep an output table
names = cell(length(whichPts),1);
car_rates = nan(length(whichPts),1);
bipolar_rates = nan(length(whichPts),1);

count = 0;
for p = whichPts
    count = count + 1;
    pt_name = pt(p).name;
    
    
    for im = 1:2
        
        fprintf('\nDoing %s\n',pt_name);
        

        %% Load spike file
        out = load([spike_folder,sprintf('%s_pc.mat',pt_name)]);
        out = out.pc;
        
        %% Skip incomplete pts
        
        % Get corresponding pt
        for j = 1:length(pt)
            if strcmp(pt(j).name,pt_name)
                break
            end
        end
        

        %% concatenate all spikes into one long thing
        % Include an extra column for the file index and block
        file_nruns = nan(length(out.file),1);
        file_spikes_per_elecs = nan(length(out.file),1);
        
        for f = 1:length(out.file)
            nspikes = 0;
            nruns = 0;
            
            
            nchs = length(out.file(f).run(1).data.montage(1).labels);
            for h = 1:length(out.file(f).run)
                nspikes = nspikes + size(out.file(f).run(h).data.montage(im).spikes,1);
                nruns = nruns + 1;
            end
            
            
            % turn to nspikes per elec
            nspikes_per_elec = nspikes/nchs;
            
            file_nruns(f) = nruns;
            file_spikes_per_elecs(f) = nspikes_per_elec;
        end
        
        %% Calculate average rate in spikes/elecs/min
        
        % weighted average across files
        avg_spike_rate = sum((file_nruns .* file_spikes_per_elecs))/sum(file_nruns)/sum(file_nruns);
        
        
        names{count} = pt_name;
        if im == 2
            car_rates(count) = avg_spike_rate;
        elseif im == 1
            bipolar_rates(count) = avg_spike_rate;
        end
        
    end
    
 
end

T = table(name,bipolar_rates,car_rates);
writetable(T,[out_folder,'spike_rates']);

end