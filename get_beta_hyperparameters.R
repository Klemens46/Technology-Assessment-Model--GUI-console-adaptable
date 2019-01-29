
# Name of script:	get_beta_hyperparameters.R
# Purpose:	A script that creates two functions to quickly calculate the Beta 
#			distribution hyper-parameters (2nd and 3rd arguments as used in 
#			the R software) from a mean and a standard deviation (SD) or a 
#			mean and a relative standard deviation (RSD).
#			This script allows to quickly parameterize a Beta 
#			distribution and if needed automatically re-parameterize it 
#			when the mean and/or the SD is changed.
# Author(s):	Klemens Wallner (who is grateful to Benoit Kudinga who 
#				spotted errors in the earlier version of this script and 
#				helped me make it more clear.)
# Note(s):		1) a = alpha, b = beta, mean = alpha / (alpha + beta)


# Function 1
# Get the Beta hyper-parameters from a mean and a SD 
Get.BHP <- function(mean, SD) {
	a <- ((1 - mean) / (SD ^ 2) - (1 / mean)) * (mean ^ 2)
	b <- a * ((1 / mean) - 1)
	return(c(a, b))
}
# For example: to draw 5 samples from a Beta distribution with a mean 
# of 0.80 and a SD of 0.1 type:
# H <- Get.BHP(0.8, 0.1)
# samples <- rbeta(5, H[1], H[2])


# Function 2
# Get the Beta hyper-parameters from a mean and a RSD
Get.BHP.RSD <- function(mean, RSD=10) {
	RSD <- RSD / 100
	SD <- mean * RSD
	a <- ((1 - mean) / (SD ^ 2) - (1 / mean)) * (mean ^ 2)
	b <- a * ((1 / mean) - 1)
	return(c(a, b))
}
# For example: to draw 5 samples from a Beta distribution with a mean 
# of 0.80 and a RSD of 7% type:
# H <- Get.BHP(0.8, 7)
# samples <- rbeta(5, H[1], H[2])
