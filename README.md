Please feel free to contact me if you have any questions or concerns.

Methods overview:

These files contain a computer simulation model. 
It's a decision analytic model and its core is a probabilistic model that calculates the standard outputs of a modern cost-effectiveness analysis. The scripts itself mention some research papers on the methods that we used. 

Further it also simulates a trial and also patient pathways after entering the hospital and after discharge at home.   
The analysis is designed to estimate the gross monetary value of a clinical trial as well as costs and benefits per patient.  This part that is estimating the gross trial value needs further work before being reliable enough to use as help in makind funding decisions.


Notes on running of these R software scripts:

The script '1_master.R' has the base case model in it and needs to be run first. 
Then comes the script '2_compute_enpvsi.R', which has the ENPVSI model in it. 
All other scripts are called from within those two initial scripts. 

If needed, the additional script 'save_other_e_results.R' saves other results of the ENPVSI model as files.
