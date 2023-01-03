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

1) run_mt_pipeline (in analyses/new_outcome/pipeline_mesial_temporal/)
 - this calls mt_patient_stitch, which in turn calls individual_run_mt
 - individual_run_mt downloads an individual edf file containing 10 minutes
 of EEG data. It first 


%}