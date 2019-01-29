
# Name of script:	get_gamma_hyperparameters.R
# Purpose:	A script that creates two functions to quickly calculate 
#			the Gamma distribution hyper-parameters shape and rate (2nd and 3rd
#			arguments of the Gamma functions as used in the R software) 
#			from a mean and a standard deviation (SD) or a mean and a 
#			relative standard deviation (RSD).
#			This script allows to quickly parameterize a Gamma 
#			distribution and if needed automatically re-parameterize it 
#			when the mean and/or the SD is changed.
# Author(s):	Klemens Wallner (who is grateful to Benoit Kudinga who helped 
#				to make this script more clear).
# Note(s):		1) a = alpha = shape, b = beta = rate , mean = shape / rate
#				2) The gamma distribution can be parameterized 
#				   using either the shape and rate (as here) or the 
#				   shape and scale parameters, where the scale = 1 / rate and 
#				   rate = 1 / scale.
#				   For example 'pgamma(1,3,2)' is equivalent to 
#				   'pgamma(1,3, scale=0.5)'.


# Function 1
# Get the Gamma hyper-parameters from a mean and a SD 
Get.GHP <- function(mean, SD) {
	b <- mean / (SD ^ 2)
	a <- mean * b
	return(c(a, b))
}
# For example: to draw 5 samples from a Gamma distribution with a mean 
# of 60 and a SD of 3 type:
# H <- Get.GHP(60, 3)
# samples <- rgamma(5, H[1], H[2])


# Function 2
# Get the Gamma hyper-parameters from a mean and a RSD in percent
# (default RSD is 10%)
Get.GHP.RSD <- function(mean, RSD=10) {
	RSD <- RSD / 100
	b <- mean / ((mean * RSD) ^ 2)
	a <- mean * b
	return(c(a, b))
}
# For example: to draw 5 samples from a Gamma distribution with a mean 
# of 60 and a RSD of 7% type:
#	H <- Get.GHP.RSD(60, 7)
#	samples <- rgamma(5, H[1], H[2])
