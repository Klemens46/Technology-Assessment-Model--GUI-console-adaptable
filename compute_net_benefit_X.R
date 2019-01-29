
# Name of script: compute_net_benefit_X.R
# Purpose: 	Calculating the monetary net benefit given the in-trial simulation results, 
#			it is a simplified version of 'compute_net_benefit.R' using 'inner.X'.
# Author(s):	Peter Hall, Klemens Wallner


# Net Benefit function ('inner.N' x "record types" x "lambda levels")
# What is recorded in 'NB': 
# The row number is the simulation (inner loop) number
# The 1st column is the individual net benefit for standard treatment
# The 2nd column is the individual net benefit for intervention
# The 3rd column shows if intervention is optimal (1) or not (0)
# The 4th column shows the maximum NB of both alternatives
# The 3rd dimension indicates the index of the lambda value used

# NB (!!): Only the 2nd column is updated here! Everything else is the same
# 		   as in the base model run.
NB.X <- NB

Compute.net.benefit.X.uncompiled <- function() {

	for (i in 1:2) {
		NB.X[, 2, i] <- inner.X[, 3] * lambda[i] - inner.X[, 4]
	}
return(NB.X)
}
