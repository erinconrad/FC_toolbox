function mats_for_nishant

HUP_id = [179  181  182  184  187  195  196  197  201  202  204  206  208  214  215  217  220];

%% File locations
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];


%% Load out file and get roc stuff
out = load([data_folder,'main_out.mat']);
out = out.out;



for i = 1:length(HUP_id)
    name = sprintf('HUP%d',HUP_id(i));
    nishant(i).name = name;
    found_it = 0;
    
    for ip = 1:length(out.all_names)
        if strcmp(out.all_names{ip},name)
            found_it = 1;
            break
        end
        
    end
    
    if ~found_it
        fprintf('\nWarning, did not find %s\n',name);
        continue;
    end
    
    % get stuff for this patient
    fc = out.all_fc{ip};
    labels = out.all_labels{ip};
    
    nishant(i).fc = fc;
    nishant(i).labels = labels;
    
end

save([plot_folder,'nishant.mat'],'nishant')

end