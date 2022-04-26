function mats_for_alfredo

%% File locations
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];


%% Load out file and get roc stuff
out = load([data_folder,'main_out.mat']);
out = out.out;
npts = length(out.all_names);


for ip = 1:npts
    
    % get stuff for this patient
    alfredo(ip).name = out.all_names{ip};
    alfredo(ip).labels = out.all_labels{ip};
    alfredo(ip).spikes = out.all_spikes{ip};
    alfredo(ip).locs = out.all_locs{ip};
    
end

save([plot_folder,'alfredo.mat'],'alfredo')

end