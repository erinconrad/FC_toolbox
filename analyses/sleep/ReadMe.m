%{

To do:
- add code to identify and track what happens to the different spike
populations over time
    - spike co-activation networks
    - identify the leader of each network, see what happens to those with
    sleep and after seizures

This contains the script to run the sleep analysis. A brief overview of
the pipeline:

0) You can add more ieeg.org patients as needed by running
create_pt_struct/add_more_pts.m
   - you must specify in the code which patients to add
   - it will find them, add their electrode locations, and then pick run
   times
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
over to borel).
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
10) To do circadian analysis, run analysis/sleep/all_pt_psd.m
11) To do analyses related to sleep, run
analysis/sleep/binary_ad_analyses.m
  - this begins by validating the alpha delta ratio against Erin's manual
  sleep/wake labels, and then selecting the optimal cut off point for the
  alpha delta ratio that best discriminates wake from sleep

%}