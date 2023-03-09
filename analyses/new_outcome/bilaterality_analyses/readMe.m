%{
To do the bilaterality prediction analysis:

Run:
>> [T,features] =  lr_mt
>> all_pred = mt_lr_loo(T,features)

mt_lr_loo (name is a typo) calls lt_mr_tree.m (in the classifiers folder)
to do the bagged ensemble classification.


%}

%% What this pulls from
%{
This uses an intermediate dataset containing many electrode-level or
electrode-to-electrode edge-level features. Steps to regenerate this
dataset:

0) the following things need to be updated in the pt structure and manual
validation file

In manual validation file:
- file start times
- SOZ

In pt structure
- outcomes
- sz times (make sure to run the get_all_annotations and then the code to
find possible seizure annotations, and manually check). This then requires
updating the manual validation file in order to add the sz_times to the
struct

1) save_edf (in save_edfs/)
- this pulls file start times from the Manual validation file (must be
downloaded from google sheets) and then downloads ieeg.org data
corresponding to a 12 hour stretch in the first night after the first full
day of recording, and it downsamples it and saves it as an edf.

2) get_sleep_stages (in save_edfs/)
- this runs the sleepSEEG pipeline on the edfs (this pipeline requires a 12
hour continuous overnight stretch of data)

3) run_mt_pipeline (in analyses/new_outcome/pipeline_mesial_temporal/)
 - this calls mt_patient_stitch, which in turn calls individual_run_mt
 - individual_run_mt downloads an individual edf file containing 10 minutes
 of EEG data. 



%}