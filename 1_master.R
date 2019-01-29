
rm(list=ls(all=TRUE))  # Remove all objects from current session memory (environment)
start.time <- proc.time()  # Record starting time of script in processor time.

# Name of script: 1_master.R [to run type: source('1_master.R')]
# Purpose of script:	Running a two alternative discrete state and time Markov model producing 
#						discounted cost & benefit outcomes with drawing of parameter values 
#						and evaluation of outputs/outcomes by means of code here and 
#						using code in other scripts (for these see other scripts). The model 
#						compares two treatments to prevent pressure ulcers (PUs) in hospital care.
# Author(s):	Peter Hall, Klemens Wallner
# Note(s): 	1)	The original R code was created by Peter Hall (then a breast cancer model).
#			2)	The alternative using standard PU prevention is 
#			called standard treatment (ST) and the alternative using standard therapy plus 
#			intermittent electric stimulation for PU prevention is called intervention (IT).
#			3) With the cycle length of 3 days and one year having about 365.25 days here one year
#			has 122 cycles (365.25/3 = 121.75 ~ 122).
#			4) Although the model has a lifetime horizon only the first 3 years 
#			are modelled using a Markov model and the remainder is modelled 
#			with an extension using linear extrapolation from t = 3 years. The 
#			extension basically takes population distribution at t = 3 years, 
#			which is 0 everywhere except in 'DF at home' and 'Dead',
#			and assumes they live for their remaining LE in 
#			those same states. This approach is much faster (~23 x) than using
#			a Markov model with a 3 day cycle length for the whole model 
#			horizon. E.g. the remaining time until they reach the age of 120, 
#			which would need to used then. Alternatively one could use another 
#			Markov model with a 1 year cycle length after the first 3 years.
#			Timing considerations are here more important, because - after our 
#			traditional CEA - we would like use to run the computationally 
#			intensive ENPVSI methods to simulate key trial outcomes.
#			5) Currently the outputs are stored in the working directory.


options(scipen=999)  # Practically disable scientific notation for numbers

# Set random number generator seed to have comparability between different 
# runs of this model and within the trial simulation model that is based on 
# this CEA model.
set.seed(1)


# ============ Set model structure ============
# (This includes to load packages and scripts, initialising variables etc.

# Load packages
library(tcltk)
#library(ggplot2)
library(MASS)
library(compiler)


######## 
# Set model structure parameters

S <- 8					# number of states 
cyclelength <- 3		# cycle length in days
T <- 366				# set time horizon in cycles (Markov model without extension)
cyclesperyear <- 122	# number of cycles per year
startage <- 52			# initial age of cohort (at t = 1 in years)
disc <- 0.015 			# discount rate for benefits and costs
K <- 2					# number of willingness-to-pay levels (lambdas)

# Note: Currently 'K' needs to be two or parts of the model might not run. Results can
#		be analysed with the package 'BCEA' or similar applications (see below).

# Set and name levels of lambdas
lambda <- c(50000, 100000)
lambda.names <- paste(as.character(lambda), sep="")
aux <- c(lambda/1000)
aux <- round(aux,3)  # Round lambda values to full thousands
lambda.names.short <- paste(as.character(aux), "k", sep="")

# Check if number of elements in 'lambda' differs from 'K'.
if (length(lambda) != K){
	stop("Error, model has stopped! Please check levels and number of lambdas.")
}

# Set if model should do a probabilistic run (1) or a 
# deterministic run (0 or everything else; for purposes of model calibration)
# and display the current setting.
# NOTE: 1) This affects the running of 'set_inputs.R' and the number of inner loops.
# 		2) Currently this just applies to transition probability parameters.
prob.run <- 1
if (prob.run == 1) {
	cat("\nThe model run is set to probabilistic. \n")
}else{
	cat("\nThe model run is set to deterministic. \n")
	inner.N <- 1
	evppi.run <- 0  # This is NOT the overall on/off switch for the EVPPI run! (see just below)
}

if (prob.run == 1) {

	# Set if model should do a EVPPI run (1) or a PSA run (0 or everything else).
	# NOTE: 1) The EVPPI run currently only applies to the parameter 'HR_PU'.
	evppi.run <- 1
	if (evppi.run == 1) {
		cat("\nThe model will do a EVPPI run.\n")
		cat("\nThe parameter of interest is 'HR_PU'.\n")

		# Ask if standard number of outer loops should be
		# used or a different number. 
		cat("\nIf you want the model to run 100 outer loops press Enter. \n")
		cat("Otherwise enter the number of loops and press Enter.")
		outer.evppi.N <- readline("\n")
		# Note the variable 'outer.N' is defined as the number of 
		# outer loops in the ENPVSI run.

		if (outer.evppi.N=="") {
			outer.evppi.N <- 100 # Define standard number of outer loops
		}else{
			outer.evppi.N <- as.numeric(outer.evppi.N)
			if (is.na(outer.evppi.N)== TRUE | outer.evppi.N < 1){
				stop("Data entering error!")
			}
		}
	}

	# Ask if standard number of simulations should be
	# used or a different number. 
	cat("\nIf you want the model to run 100 simulations press Enter. \n")
	cat("Otherwise enter the number of simulations and press Enter.")
	inner.N <- readline("\n")

	if (inner.N=="") {
		inner.N <- 100 # Define standard number of simulations for normal run
	}else{
		inner.N <- as.numeric(inner.N)
			if (is.na(inner.N)== TRUE | inner.N < 1){
				stop("Data entering error!")
			}
	}
}
########


########
# Load scripts and compile the model script
# Note: the script 'compute_net_benefit.R' needs to be loaded later on.

# Load Beta hyper-parameter functions from script 
#(used in 'set_inputs.R' below)
source('get_beta_hyperparameters.R')

# Load Gamma hyper-parameter functions from script 
#(used in 'set_inputs.R' below)
source('get_gamma_hyperparameters.R')

# Load 'set inputs' script i.e. define non-structural inputs 
#(e.g. draw sample values for every simulation i.e. inner loop)
source('set_inputs.R')
# Note: if 'evppi.run == 1' then each outer loop will 
# use the same inner loop draws.


# Load single model iteration function 'Run.model', which 
# returns 'c(qaly1, cost1, qaly2, cost2)'
source('run_model.R')

# Compile 'model' function for faster running
Run.model <- cmpfun(Run.model.uncompiled)
# Note: For the other functions it does not 
# pay off to compile them for 'inner.N' <= 10,000 because the compiling 
# takes longer than the reduction in running time.
########


########
#Initialise variables and create progress bar
# Note: this process differs between a normal and a EVPPI run.


if (evppi.run != 1){

	# Variable to record cost-effectiveness results of inner loops (simulations)
	# Note: these are currently in this order: QALY.s, COST.s, QALY.i, COST.i, LY.s, LY.i
	inner <- array(c(NA, NA), c(inner.N, 6))

	# Check variables record for model calibration
	# These can vary - see return function of 'run_model.R'
	ch <- array(NA, c(inner.N, 12))

	####
	# Create progress bar and display that progress bar is for simulations part only
	cat("\nProgress bar for the simulation part of the model:\n")
	pb <- txtProgressBar(min = 0, max = inner.N, style = 3)
	####

	# Vector that shows for each simulation if there was an 
	# error with the row sums (i.e. is the sum of each TP row = 1?) yes (1) or 
	# no (0) .
	# Note: 1) If any of the TP entries is smaller than 0 the whole 
	#		model is stopped.
	#		2) If 1) is not the case and the row sum is > 1 
	#		then the row sum error is fixed right after detection by normalising
	#		the TP entries (for details see script 'run_model.R').
	rowSumsError_dummy <- rep(NA, inner.N)

}else{
	# Variable to record cost-effectiveness results of outer loops (samples)
	outer.loops <- array(c(NA, NA), c(outer.evppi.N, 6))

	####
	# Create progress bar and display that progress bar is for simulations part only
	cat("\nProgress bar for the simulation part of the model:\n")
	pb <- txtProgressBar(min = 0, max = outer.evppi.N, style = 3)
	####

	# Note: This might work differently in the EVPPI run.
	rowSumsError_dummy <- rep(NA, outer.evppi.N)

	# Variable to record the CE results of the inner EVPPI loops (simulations)
	inner.evppi <- array(c(NA, NA), c(inner.N, 6)) 
}


# temporary record of all model returns
MR <- rep(NA, 18)
########

# FYI - currently the model return from 'run_model.R' is: 
#model_return <- c(QALY.total.s, COST.total.s, 
#QALY.total.i, COST.total.i, LY.s, LY.i, cv, COST.s.rem, COST.i.rem, QALY.s.rem, QALY.i.rem)


# ============ Run simulation model ============
if (evppi.run != 1){
	# Normal PSA model run (also used for "deterministic" run).
	for (i in 1:inner.N) {
		# Record model return (MR) (return of the model function i.e. recorded results)
		MR <- Run.model(i)
		inner[i, ] <- MR[1:6]
		ch[i, ] <- MR[7:18]
		# Update progress bar
		setTxtProgressBar(pb, i)
	}
}else{
	# EVPPI model run
	HR_PU_outer <- HR_PU
	for (o in 1:outer.evppi.N) {
		# Define parameter of interest for this outer loop
		HR_PU <- rep(HR_PU_outer[o], inner.N)

		for (i in 1:inner.N){
		# Record model return (MR) (return of the model function i.e. recorded results)
		MR <- Run.model(i)
		inner.evppi[i, ] <- MR[1:6]
		}

		# Record expected values of this outer loop.
		outer.loops[o,] <- colMeans(inner.evppi)

		# Update progress bar
		setTxtProgressBar(pb, o)
		
	}
	HR_PU <- HR_PU_outer
}


if (evppi.run == 1){
	# In order to make code below work correctly for an EVPPI run we need to 
	# re-define 'inner' and 'inner.N'.
	inner <- outer.loops
	inner.N.save <- inner.N
	inner.N <- outer.evppi.N
}


# ============ Calculate net monetary benefit (NB) for simulations ============

# Set expanded features in 'compute_net_benefit.R' to on (1) or off (0).
# Their results are needed for the ENPVSI calculation.
X <- 1

# Load net benefit function
# Note: This loading of the net benefit function needs to be at this
# exact code location. 
source("compute_net_benefit.R")

NB <- array(NA, c(inner.N, 4, 2))
NB <- Compute.net.benefit.uncompiled()

# ============ EVPPI run only: redefine 'inner' and 'inner.N' ============


# ============ Calculate model outputs ============

c.Sta <- mean(inner[, 2])  # mean costs standard treatment
c.Int <- mean(inner[, 4])  # mean costs intervention
e.Sta <- mean(inner[, 1])  # mean QALYs standard treatment
e.Int <- mean(inner[, 3])  # mean QALYs intervention

e.LY.Sta <- mean(inner[, 5])  # mean LYs standard treatment
e.LY.Int <- mean(inner[, 6])  # mean LYs intervention


c.inc.Int <- c.Int - c.Sta  # mean incremental costs of intervention
e.inc.Int <- e.Int - e.Sta  # mean incremental QALYs of intervention
e.LY.inc.Int <- e.LY.Int - e.LY.Sta  # mean incremental LYs of intervention

# Incremental cost-effectiveness ratio (QALY)
ICER <- c.inc.Int / e.inc.Int  

########
# Expected NB with current information

NB.current <- NB.1  <- NB.2 <- rep(NA, K)

for (k in 1:K) {

# Maximum expected NB for each lambda level (one value for each level)
	NB.current[k] <- max(mean(NB[, 1, k]), mean(NB[, 2, k]))
	
# Expected NB for each alternative
	NB.1[k] <- mean(NB[, 1, k])
	NB.2[k] <- mean(NB[, 2, k])
}

# Create matrix with knowledge about expected NBs with current information:
# 1. col: NBs of strategy 1; 2. col: NBs of strategy 2, 
# 3. col: optimal strategy at this lambda; 4. col: NB of optimal strategy
# (This matrix will later be used as the current information reference).
NBs.current <- matrix(NA, nrow=K, ncol=4)
NBs.current[, 1] <- NB.1
NBs.current[, 2] <- NB.2
NBs.current[, 4] <- NB.current
for (k in 1:K ) {
	NBs.current[k, 3] <- if (NBs.current[k, 2] >= NBs.current[k, 1]){
							2
						}else{
							1
						}
}
colnames(NBs.current) <- c("NB treatment 1", "NB treatment 2", "Optimal", "Maximum NB")
rownames(NBs.current) <- c(lambda.names)
write(NBs.current, "NBs.current.data")
########


########
# Expected value of perfect information (EVPI) calculation
EVPI  <- array(c(lambda, rep(NA, 2)), c(K, 2))
for (k in 1:K) {	
	EVPI[k, 2] <- mean(NB[, 4, k]) - NB.current[k]
}
########


# ============ Display model outputs ============

# Mean outputs for each treatment alternative
ce.mean.outputs <- matrix(c(e.Sta, c.Sta, e.Int, c.Int), nrow=2, ncol=2)
colnames(ce.mean.outputs) <- c("Standard treatment (1)", "Intervention (2)")
rownames(ce.mean.outputs) <- c("QALYs", "Cost")
cat("\n\nMean results for each of the two alternatives")
cat("\n (the column number indicates the alternative) \n")
print(ce.mean.outputs)

# Currently not always working correctly:
effect <- matrix(c(inner[, 1], inner[, 3]), nrow=inner.N, ncol=2)
cost   <- matrix(c(inner[, 2], inner[, 4]), nrow=inner.N, ncol=2)
cat("\n\n The two objects ('e' and 'c') can be created (see below) and then\n")
cat(" analysed with the package 'BCEA', once installed. \n")
cat(" Note: The package 'BCEA' uses some unusual object names, \n please be aware of this. \n\n")
cat("Creation of 'e' and 'c':\n e <- effect\n")
cat(" c <- cost\n")


# Categorize the simulations 
# into the 5 categories (= 4 CE-plane quadrants plus origin) compared to standard treatment:
# - more costly and more effective (north east - NE)
# - more costly and equally or less effective, or 
#	equally costly and less effective (north west - NW, i.e. dominated)
# - less costly and less effective (south west - SW)
# - less costly and equally or more effective, or 
#	equally costly and more effective (south east - SE, i.e. dominating)
# - equally effective and equally costly (origin - OR)
# Notes:	1)This section is currently only useful for 2 alternatives and always
# 			uses alternative 1 as reference case.
#			2) The borders of the SE quadrant logically belong to the SE quadrant (dominating) and
#			the borders of the NW quadrant logically belong to the NW quadrant (dominated).

# Create a matrix that shows the categorisation above.
ce.plane.cat <- matrix(0, nrow=inner.N, ncol=2)
ce.plane.cat[, 1] <- "NA"  # Standard treatment is the reference
for (s in 1:inner.N) {
	ce.plane.cat[s, 2] <-	if ((inner[s, 3] > inner[s, 1]) & (inner[s, 4] > inner[s, 2])) {
								"ne"
							}else{
								if (((inner[s, 3] >= inner[s, 1]) & (inner[s, 4] < inner[s, 2])) || ((inner[s, 3] > inner[s, 1]) & (inner[s, 4] == inner[s, 2]))){
									"se"
								}else{
									if (((inner[s, 3] <= inner[s, 1]) & (inner[s, 4] > inner[s, 2])) || ((inner[s, 3] < inner[s, 1]) & (inner[s, 4] == inner[s, 2]))) {
										"nw"
									}else{
										if((inner[s, 3] < inner[s, 1]) & (inner[s, 4] < inner[s, 2])) {
											"sw"
										}else{
											if((inner[s, 3] == inner[s, 1]) & (inner[s, 4] == inner[s, 2])) {
											"or"
											}
										}
									}
								}
							}
}

# Calculate and display the number of simulations in each category for alternative 2.
# Not used at the moment:
sum.ne <- sum(ce.plane.cat[, 2] == "ne")
sum.nw <- sum(ce.plane.cat[, 2] == "nw")
sum.se <- sum(ce.plane.cat[, 2] == "se")
sum.sw <- sum(ce.plane.cat[, 2] == "sw")
sum.or <- sum(ce.plane.cat[, 2] == "or")
cat("\nThe", c(inner.N), "simulations of treatment 2 fall into")
cat("\nthe following 5 categories of the CE-plane (using QALYs):\n\n")
cat(" - More costly and more effective (NE):", c(sum.ne), "\n")
cat(" - More costly and equally or less effective, or \n   equally costly and less effective (NW + borders, i.e. dominated): ", c(sum.nw), "\n")
cat(" - Less costly and less effective (SW): ", c(sum.sw), "\n")
cat(" - Less costly and equally or more effective, or \n   equally costly and more effective (SE + borders, i.e. dominating): ", c(sum.se), "\n")
cat(" - Equally effective and equally costly (origin - OR): ", c(sum.or), "\n")


# Include code to save the CEA summary in a separate file and display it:

# Create matrix that contains a CEA summary.
# Note: Best CE decision at different WTP levels: 0 = is not best, 1 = is best option.
CEA.summary <- matrix(c(ce.mean.outputs), nrow=2, ncol=2)
aux.matrix <- matrix(c( e.LY.Sta, e.LY.Int,
						- e.inc.Int, e.inc.Int, 
						- c.inc.Int, c.inc.Int, 
						- e.LY.inc.Int, e.LY.inc.Int, 
						ICER, ICER, 
						NBs.current[1, 1], NBs.current[1, 2], 
					   NBs.current[2, 1], NBs.current[2, 2], 
					   NBs.current[1, 1] - NBs.current[1, 2], # Part 1 of matrix row 6
					   NBs.current[1, 2] - NBs.current[1, 1], # Part 2 of matrix row 6
					   NBs.current[2, 1] - NBs.current[2, 2], # Part 1 of matrix row 7
					   NBs.current[2, 2] - NBs.current[2, 1], # Part 2 of matrix row 7
					   if (NBs.current[1, 3] == 1) {
							1
                       } else {
							0
                       }
					   , 
					   if (NBs.current[1, 3] == 2) {
							1
                       } else {
							0
                       }
					   , 
                       if (NBs.current[2, 3] == 1) {
							1
                       } else {
							0
                       }
					   , 
                       if (NBs.current[2, 3] == 2) {
							1
                       } else {
							0
                       }
					   , 
					   c(sum.sw * 100 / inner.N), c(sum.ne * 100 / inner.N), 
					   c(sum.se * 100 / inner.N), c(sum.nw * 100 / inner.N), 
                       c(sum.ne * 100 / inner.N), c(sum.sw * 100 / inner.N), 
                       c(sum.nw * 100 / inner.N), c(sum.se * 100 / inner.N), 
		               c(sum.or * 100 / inner.N), c(sum.or * 100 / inner.N), 
					   EVPI[1, 2], EVPI[1, 2], 
					   EVPI[2, 2], EVPI[2, 2]
					   ), nrow=18, ncol=2, byrow=T)
CEA.summary <- rbind(CEA.summary, aux.matrix)
colnames(CEA.summary) <- c("Treatment 1", "Treatment 2")
rownames(CEA.summary) <- c("Mean QALYs", "Mean cost", "Mean LYs (un-discounted.)",
							 "Incr. QALYs", "Incr. cost", "Incr. LYs (un-discounted.)", 
							 "ICERs(!)", 
							 paste("Indiv. NMB at WTP=", c(lambda.names.short[1]), sep=""), 
						     paste("Indiv. NMB at WTP=", c(lambda.names.short[2]), sep=""), 
							 paste("Incr. NMB at WTP=", c(lambda.names.short[1]), sep=""), 
						     paste("Incr. NMB at WTP=", c(lambda.names.short[2]), sep=""), 
							 paste("Best decision at WTP=", c(lambda.names.short[1]), sep=""), 
							 paste("Best decision at WTP=", c(lambda.names.short[2]), sep=""),
						     "% More costly and more effective (NE):", 
						     "% Dominated (NW + borders):", 
						     "% Less costly and less effective (SW):", 
						     "% Dominating (SE + borders):", 
							 "% Equally costly and equally effective (OR):", 
						     paste("EVPI per patient at WTP=", c(lambda.names.short[1]), sep=""),
						     paste("EVPI per patient at WTP=", c(lambda.names.short[2]), sep=""))

# Saving the matrix 'CEA.summary' in a text file. One can use the argument 
# 'append=T' meaning that if there is already a file with the same name, 
# the matrix will be saved below any saved content (this way can however 
# cause a warning message).
write.table(CEA.summary, "CEA-summary.txt")

# Saving the matrix 'CEA.summary' in a CSV file that can be easily opened and edited in MS Excel etc.
write.csv(CEA.summary, file=paste("CEA_summary.csv", sep=""), row.names=FALSE)


# Display the 'CEA.summary' matrix and suggest to show 3 digits and 
# leave at least a 4 spaces gap between columns (default=1).
cat(" \n", "The following summary of this Cost-Effectiveness Analysis was")
cat(" saved in \n the files 'CEA_summary.csv' and 'CEA-summary.txt': \n \n")
print(CEA.summary, digits=5, print.gap=4)
# Note: R may display much more significant digits if the results are 
# presented in a table and one of the elements of the column needs a lot of 0 before its significant digits. 
# In that case R adjusts the all values in that column to the same number of digits (significant or not) so that
# the corresponding digits and dots of each element are exactly below each other. For instance:
# 25.0000
#  0.0012

options(digits=10)

########
# Show number of loops/simulations and if it's not a EVPPI-run of the model also show the 
# outputs needed for calibration and report on TP errors 
if (evppi.run != 1){
	cat("\n\nNumber of simulations:", c(inner.N), "\n")

	###
	# Show outputs needed for calibration (currently still for an model)
	# Note: currently still from an older model!
	cat("\n\n Mean duration in Hospital (in days, 1, 2, 3, 4; standard care)\n", c(mean(ch[, 1])))
	
	cat("\n\n Mean duration of being with 1+ PU (in days per real PU patient; standard care)\n", c(mean(ch[, 2])/(sum(mean(ch[, 4]),mean(ch[, 6])))))
	
	
	#cat("\n\n [No current description of most calibration values entered.]", c(mean(ch[, 2])))
	#cat("\n\n ", c(mean(ch[, 3])))  # Currently not used
	cat("\n\n Incidence of >=1 PU (standard care) ", c(mean(ch[, 4])))
	cat("\n\n Incidence of >=1 PU (intervention) ", c(mean(ch[, 5])))

	cat("\n\n Incidence of = 2 PU (standard care)", c(mean(ch[, 6])))
	#cat("\n\n Ratio of patients that develop
	# a 2nd PU (only from PU2 or A0; HSF arm)", c(mean(ch[, 7])))
	#cat("\n\n Ratio of patients that develop  # Currently not used
	#a 2nd PU (only from PU2 or A0; APM arm)", c(mean(ch[, 8])))

	# Put row sum error information in the correct object
	rowSumsError_dummy <- ch[, 8]

	# Error report for transition probability errors
	if(sum(rowSumsError_dummy) == 0){
		cat("\n\nNo errors in any of the transition probabilities.\n")
	}else{
		cat("\n\nNumber of row sum errors in the transition probabilities:\n")
		cat(c(sum(rowSumsError_dummy)))
	}
	###

}else{
	cat("\n\nNumber of outer loops: ", c(outer.evppi.N))
	cat("\n\nNumber of inner loops: ", c(inner.N.save), "\n")
}
########


# Display Net benefit results
cat("\n Net benefit results:\n")
print(NBs.current)


###
# Save workspace using unique name (to reuse objects later)
uname <- format(Sys.time(), "%Y_%b_%d-%H_%M_%S")
save(list = ls(all = TRUE), file = paste("post_1_master-",
	 as.character(inner.N), "-", as.character(uname), ".RData", sep=""))
cat("\n Workspace was stored after running the model \n")
cat("(in file: 'post_1_master-", c(inner.N), "-",c(uname), ".RData')\n", sep="")
###


########
## Display a scatter plot of the cost-effectiveness plane
#e.inc.Int.ps <- inner[,3] - inner[,1]
#c.inc.Int.ps <- inner[,4] - inner[,2]
#
#xlimit <- max(abs(e.inc.Int.ps))*1.1
#ylimit <- max(abs(c.inc.Int.ps))*1.1
#
#plot(e.inc.Int.ps, c.inc.Int.ps, xlab="Incremental Benefit (QALY)",
# ylab="Incremental Cost ($)", xlim=c(-xlimit,xlimit), ylim=c(-ylimit,ylimit),
# main="Cost-Effectiveness Plane", pch=20, cex.main=1.5, 
# frame.plot=T, col="blue")
#
#abline(v=0, col="black")
#abline(h=0, col="black")
#########


# ============ Final timing part ============  

# Suggest the R software to display 2 digits.
options(digits=2)

# Calculate running time of script by subtracting processor time at the start
# from current processor time.
duration <- matrix(c(start.time - proc.time()))

# Display the time it took to run the script 
cat("\n \n Running this script took ", c( - duration[3, ]), "seconds ")
cat("(~", c(( - duration[3, ]) / 60), "minutes). \n")

options(digits=5)  # Reset digits display to 5

# Scatter plot on CE plane :-)
#
e.inc.Int.ps <- inner[,3] - inner[,1]
c.inc.Int.ps <- inner[,4] - inner[,2]


plot(e.inc.Int.ps, c.inc.Int.ps, xlab="Incremental Benefit (QALY)",
 ylab="Incremental Cost ($)", 
 main="Cost-Effectiveness Plane", pch=20, cex.main=1.5, 
 frame.plot=T, col="blue")
#
#abline(v=mean(e.inc.Int.ps, na.rm=TRUE), col="orange")
#abline(h=mean(c.inc.Int.ps, na.rm=TRUE), col="red")
#
abline(v=0, col="black")
abline(h=0, col="black")

#
# To comment a whole block out: in the selected block, replace \n with \n# and then 
# put a # in the first row of the block

