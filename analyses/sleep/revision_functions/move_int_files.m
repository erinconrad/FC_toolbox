function move_int_files

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
old_epilepsia_folder = [results_folder,'analysis/backup_intermediate_Feb26_good_spikes/'];
new_all_folder = [results_folder,'analysis/intermediate/'];
new_epilepsia_folder = [results_folder,'analysis/intermediate_epilepsia_revision/'];

%% Loop over all intermediate
listing = dir([new_all_folder,'*mat']);
for l = 1:length(listing)
    fname = listing(l).name;

    % see if this patient exists in the old epilepsia folder
    if exist([old_epilepsia_folder,fname],'file') ~= 0
        
        % Copy the file into the new epilepsia folder
        copyfile([new_all_folder,fname],[new_epilepsia_folder,fname]);
    end

end


end