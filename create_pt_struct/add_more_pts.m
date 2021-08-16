% Specify which hup nums to add
first_and_last_ieeg_nums = [99 99];

% Get pt info from ieeg.org and add to pt struct
fprintf('\nGetting ieeg.org info\n');
find_intracranial_pts(first_and_last_ieeg_nums);

% electrode locs
fprintf('\nGetting electrode locs\n');
add_elec_locs

% define run times
fprintf('\nGetting run times\n');
define_run_times