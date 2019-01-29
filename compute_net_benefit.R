
# Name of script: compute_net_benefit.R
# Purpose: Calculating the individual net monetary benefit (NB) given the simulation results.
# Author(s):	Peter Hall, Klemens Wallner

# Net Benefit function ('inner.N' x "record types" x "lambda levels")
# What is recorded in 'NB': 
# The row number is the simulation (inner loop) number
# The 1st column is the individual net benefit for standard treatment
# The 2nd column is the individual net benefit for intervention
# The 3rd column shows if intervention is optimal (1) or not (0)
# The 4th column shows the maximum NB of both alternatives
# The 3rd dimension indicates the index of the lambda value used

Compute.net.benefit.uncompiled <- function() {

# X is the variable that shows if the expanded features are switched on (1) or off (0).
# Switching them off reduces running time, when expanded features are not needed.
if (X == 1) {
	for (i in 1:2) {
		NB[, 1, i] <- inner[, 1] * lambda[i] - inner[, 2]
		NB[, 2, i] <- inner[, 3] * lambda[i] - inner[, 4]

			for (j in 1:inner.N) {
				NB[j, 3, i] <- if (NB[j, 2, i] > NB[j, 1, i]) 1 else 0
				NB[j, 4, i] <- max(NB[j, 1, i], NB[j, 2, i])
			}
	}
} else {
	for (i in 1:2) {
		NB[, 1, i] <- inner[, 1] * lambda[i] - inner[, 2]
		NB[, 2, i] <- inner[, 3] * lambda[i] - inner[, 4]
	}
}
return(NB)
}

