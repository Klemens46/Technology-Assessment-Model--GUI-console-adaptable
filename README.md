Please feel free to contact me if you have any questions or concerns.

Methods overview:

These files contain a two-part computer simulation model.
It's a decision analytic model and its core is a probabilistic model. This calculates the standard outputs of a modern cost-effectiveness analysis (part 1). It also aims to estimate the value of running a trial comparing two different technologies (part 2). The scripts itself mention research papers on the methods that are used. The core model is a state transition model. 

The second part simulates using one technology alongside a trial of new and current technologies. 
The analysis is designed to estimate the gross monetary value of a clinical trial and costs and benefits per patient. This analysis includes technology life time, sample size, and other variables. The second part needs further work before being reliable enough to use as help in making real-life decisions.


Notes on running the R software scripts:

The script '1_master.R' has the base case model in it and needs to be run first. 
Then comes the script '2_compute_enpvsi.R', which has the trial simulation model in it. 
All other scripts are called from within those two initial scripts. 

If needed, the additional script 'save_other_e_results.R' saves other results of the trial model as files.

I tried to comment the code as good as possible. 
Again, please feel free to contact me if you have any questions or concerns.
