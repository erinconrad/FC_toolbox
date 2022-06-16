%{
This folder contains the scripts to perform the analyses from the
Intracranial EEG Network Biases paper. 

The code for this analysis can be found at:
https://github.com/erinconrad/FC_toolbox

%}

%% To re-generate plots and tables and re-run statistical tests
%{
To run the code to regenerate the results plots and tables, perform the
following steps:
1) Download the codebase from https://github.com/erinconrad/FC_toolbox/
2) Download the intermediate datasets located in Erin Conrad's Penn Box:
https://upenn.box.com/s/xsins7ud83e4y93w3yf5uzm37ey40ir4
3) Add a file named fc_toolbox_locs.m to your Matlab path. This points to
various other paths. Here is an example of what it would look like:
    function locations = fc_toolbox_locs
        locations.paper_plot_folder = [path to where you will output results from the
            analysis]
        locations.paper_data_folder = [path to where you put the intermediate datasets]
        locations.script_folder = [path to this code base]
    end
4) Navigate to the [locations.script_folder,'analyses/outcome_prediction/plots/data/'] folder
5) Run the following command:
>>  main_ieeg_biases
This function takes an optional argument equal to the number of
testing/training splits for the SOZ classifier. If left blank then this
will default to 1,000 (the number used in the paper). It will take several
hours to run this with 1,000 splits. To run it in a few minutes, try 10-20
splits.

The takes intermediate datasets that contain functional networks, spike
rate data, at atlas parcellations for each patient. It runs the main 
analyses and generates figures and the results section text.
%}

%% Info
%{
Erin Conrad
University of Pennsylvania
2022
%}
