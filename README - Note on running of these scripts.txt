
The script '1_master.R' has the base case model in it and needs to be run first. 

Then comes the script '2_compute_enpvsi.R', which has the ENPVSI model in it. All other scripts are called from within those two initial scripts. 

If needed the additional script 'save_other_e_results.R' saves other results of the ENPVSI model as files.

