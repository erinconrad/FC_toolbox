%{
To do the bilaterality prediction analysis:

Run:
>> [T,features] =  lr_mt
>> all_pred = mt_lr_loo(T,features)

mt_lr_loo (name is a typo) calls lt_mr_tree.m (in the classifiers folder)
to do the bagged ensemble classification.


%}