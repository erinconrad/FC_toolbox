function add_more_pts(which_step,which_pts)

switch which_step
    
    case 'initial'
        
        
        fprintf(['\nEnsure that you have done the following things:\n',...
            '- added file start times to manual validation\n',...
            '- Added seizure times and SOZ to Manual Validation file\n',...
            '- saved Manual Validation to scripts/spike_detector\n']);
        pause
        
        % Specify which hup nums to add
        first_and_last_ieeg_nums = which_pts;

        % Get pt info from ieeg.org and add to pt struct
        fprintf('\nGetting ieeg.org info\n');
        find_intracranial_pts(first_and_last_ieeg_nums);

        % electrode locs
        fprintf('\nGetting electrode locs\n');
        add_elec_locs

        % define run times
        fprintf('\nGetting run times\n');
        define_run_times

        % add file start times
        add_file_start_times

        % add demographics
        add_demographics

        % mine annotations for seizures
        mine_annotations_for_szs
        
        % Add seizure times
        add_sz_times_multiple_sources
        
    case 'post_manual_addition'
        
        fprintf(['\nEnsure that you have done the following things:\n',...
            '- Done spike and FC calculations: do_run/long_run and do_run/run_for_coherence)\n',...
            '- Plotted spike detections (scripts/spike_detector/plot_example_detections)\n',...
            '- Updated Manual validation file with PPVs for spike detections\n',...
            '- saved Manual Validation to scripts/spike_detector\n',...
            ]);
        
        pause
        
        fprintf(['\nOnce all of this is done, you should be ready to get intermediate dataset '...
            '(run scripts/analysis/prep_out_file/intermediate_data\n']);
        
        

end