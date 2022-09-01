%% Summary

%{
This contains the script to run the spikes and sleep analysis. 

The code for this analysis can be found at:
https://github.com/erinconrad/FC_toolbox

%}

%% To re-generate plots and tables and re-run statistical tests
%{
To run the code to regenerate the results plots and tables, perform the
following steps:

1) Download the codebase from https://github.com/erinconrad/FC_toolbox/
2) Download CircStat, which the code calls:
https://www.jstatsoft.org/article/view/v031i10
3) Add a file named fc_toolbox_locs.m to your Matlab path. This points to
various other paths. Here is an example of what it would look like:

    function locations = fc_toolbox_locs

        locations.main_folder = [path to where you will output results from the
            analysis]
        locations.script_folder = [path to this code base]

    end
4) Navigate to the analyses/sleep/Epilepsia_plots_and_tables/ folder
5) Run the following command:
>> epilepsia_plots_and_tables

This command take an intermediate dataset called out.mat which contains
summary level information on spike rates for each patient and will execute
the code to generate plots, tables, and perform statistical testing. It
does not regenerate the Supplemental Table 1 because this requires using an
alternate out.mat intermediate dataset calculated using different
peri-ictal window definitions.

%}

%% To re-run the full data analysis from scratch
%{
This will take days-to-weeks to run the spike detector on all patients'
data.

To re-run the full pipeline:

0) You can add more ieeg.org patients as needed by running
create_pt_struct/add_more_pts.m
   - you must specify in the code which patients to add
   - it will find them, add their electrode locations, pick run
   times, and add file start times and demographics
   - note that the latter two will require digging into the chart and so
   require tables that must be added manually
1) detect spikes, do network and alpha delta ratio calculations by running
do_run/long_run.m
   - specify which patients to run
2) Generate plots of 50 random spikes for each patient in
spike_detector/plot_example_detections.m
3) Update the google sheet Manual validation.xlsx with the PPV. This will
automatically update a column in this sheet containing the HUP IDs of the
good spike detection patients
4) Download this updated spreadsheet and put it in the SCRIPT folder
spike_detector/ (it's in a script folder so that it will automatically sync
over to borel, the server Erin used to run the analysis).
5) Finish adding seizure times for patients with missing data in the google
spreadsheet seizure times.xlsx. Download this into the DATA folder
data/sz_times/
6) Update seizure times in pt.mat struct by running
seizure_info/seizure_times/add_sz_times_multiple_sources.m
7) Make sure to copy the data struct pt.mat over to Borel. It will look for
seizures in this in future steps.
8) Generate intermediate files by running
analyses/prep_out_file/intermediate_data.m. This stitches spikes and
networks from different files together, gets seizure times, and removes
spikes in seizures.
9) (already done) Generate manual labels for sleep and wake
   - analyses/sleep/get_sleep_labels/plot_scalp.m: this finds patientd with
   scalp eeg data and picks a bunch of random minutes and saves the eeg
   data for plotting
   - analyses/sleep/get_sleep_labels/scalp_gui.m: this takes the saved EEG
   data above and plots the data. Erin went through and manually labeled
   the sleep state.
10) To do analyses related to sleep, run
analysis/sleep/run_sleep_analyses.m
  - this begins by validating the alpha delta ratio against Erin's manual
  sleep/wake labels, and then selecting the optimal cut off point for the
  alpha delta ratio that best discriminates wake from sleep
  - It then does various things looking at spikes and seizures related to
  sleep
  - This generates an out.mat file, which you should move to the
  scripts/analuses/sleep/data folder
11) Run analysis/sleep/plots and tables/sleep_plots_and_tables.m
  - This will take the out.mat file and run scripts to generate plots and
  do statistical tests

%}

%% Info
%{
Erin Conrad
University of Pennsylvania
2021
%}