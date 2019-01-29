start.time <- proc.time()  # Record starting time of script in processor time.

# Name of script:	2_compute_enpvsi.R [to run type: source('2_compute_enpvsi.R')]
# Purpose: 	To calculate the expected net present value of sample information (ENPVSI) calculation 
#			for the PU model (in parent script '1_master.R') for different trial sizes.
# Author(s):	Peter Hall, Klemens Wallner
# Note(s):	1) The script '1_master.R' needs to before running this script, which is
#			based on those results and objects created.
#			2) The steps from the general Monte Carlo sampling algorithm for calculation of 
#			population ENPVSI are denoted as A1. to A3., and refer to the (adjusted) algorithm in
#			the appendix of the article by Hall et al. (2011):
#			Hall, P. S., Edlin, R., Kharroubi, S., Gregory, W., & McCabe, C. (2012). Expected net 
#			present value of sample information: from burden to investment. 
#			Medical Decision Making, 32(3), E11–21. doi:10.1177/0272989X12443010.
#			3) Ades et al. (2004) refers to:
#			Ades, A. E., Lu, G., & Claxton, K. (2004). Expected value of sample information 
#			calculations in medical decision modeling. Medical Decision Making, 24(2), 
#			207–27. doi:10.1177/0272989X04263162.
#			4) Currently the outputs are stored in the working directory.


# ENPVSI2-4 = posterior net benefit for all but data only for within trial and for trial 
# simulation - correct recruitment vector in trial sim for multiple sample sizes

print(date()) # Print current time and date to record start of script.


# ============ Load scripts and set model structure ============

# Load single model iteration function 'Run.model.intervention.arm.uncompiled' for 
# intervention arm, which returns 'c(qaly2, cost2)' and will be used to simulate the
# intervention posterior distribution results for the not-in-trial group.
# Note: The standard treatment arm is not simulated again. In that arm nothing 
# changes (because it is independent of 'HR_PU' and otherwise the same distribution draws are 
# used in the EVSI analysis), therefore its results are later on set to be equal 
# to the base case run results (derived from '1_master.R').
source("run_model_intervention_arm.R")

# Set expanded features in 'compute_net_benefit.R' to on (1) or off (0) when loading.
# Here 0 because they are not directly needed for ENPVSI calculation.
X <- 0

# Reload net benefit function for the not-in-trial group (uses 'inner')
# Note: Needed because 'X' is now set to 0. This could (!) make the model run quicker.
source("compute_net_benefit.R")


# Load trial time estimation function
source("estimate_trial_time.R")

# Load single model iteration function 'Run.model.intervention.arm.X.uncompiled(i)' for 
# intervention arm, which returns 'c(qaly2, cost2)' and will be used to simulate the 
# intervention sample X based results 
# for the in trial group.
source("run_model_intervention_arm_X.R")

# Load net benefit function for the in-trial group (uses 'inner.X') 
source("compute_net_benefit_X.R")

# Loaded functions are later compiled for faster running.
# But this is done directly before the EVSI specific code because some variables are not yet
# created/initialized, which can cause problems with compiling.

########
# Define ENPVSI-specific variables of the model structure

###
# Setting the main approach:
# During the trial is the NB maximising alternative used for not-trial population?
# Note: OWR means Only-With-Research, OIR means Only-In-Research

#Initialising 'OIR' and 'OWR'
OIR <- "OIR"
OWR <- "OWR"
# Defining approach 
approach <- OIR
###


# It is important to note that 'inner.N' and other variables used but not 
# defined here, are defined in '1_master.R' or the scripts it called. This 
# includes the base case and current information results.

# Ask if standard number of simulations should be
# used or a different number.
cat("\nThe model will run", c(inner.N), "inner loops for every outer loop.\n")
cat("\nIf you want the model to run 1 outer loop  press Enter.\n")
cat("Otherwise enter the number of outer loops and press Enter.")
outer.N <- readline("\n")
if (outer.N=="") {
	outer.N <- 1  # Define number of outer loops for normal run  
}else{
	outer.N <- as.numeric(outer.N)
		if (is.na(outer.N)== TRUE | outer.N < 1){
			stop("Data entering error!")
		}
}

# Make sure that the number of outer loops is not greater than the number of inner
# loops because that would not work currently (but can be implemented).
outer.N <- min(outer.N, inner.N)

# Number of sample numbers (different trial sizes)
N <- 2

# Defining number of patients in each arm of the simulated trial.
# Note: accuracy of the model increases the further away from 1 the sample numbers are. 
N1 <- c(100, 200)	# number of patients in standard arm 
N2 <- c(100, 200)	# number of patients in intervention arm 




# Useful lifetime of invention (assumption)
#
# Ask if default figure should be
# used or a different integer.
cat("Useful lifetime of invention (assumption)")
cat("\nIf you want the model to use a value of 200, press Enter.\n")
cat("Otherwise enter the number and press Enter.")
life <- readline("\n")
if (life=="") {
	life <- 200  # Define default figure (as integer) 
}else{
	life <- as.numeric(life)
		if (is.na(life)== TRUE | life < 1){
			stop("Data entering error!")
		}
}

###
# Define relevant patient population (IP)
#
# Patient population: Number of patients in Alberta per year (IP) that could potentially benefit from new 
# information generated in the trial - cohort that could use preventative treatment.
#
# Ask if standard population number should be
# used or a different number.
IP <- NA  # initialisation
cat("\nIf you want the model to use a relevant population of 50000 press Enter.\n")
cat("Otherwise enter the number and press Enter.")
IP <- readline("\n")
if (IP=="") {
	IP <- 50000  # Define default population number
}else{
	IP <- as.numeric(IP)
		if (is.na(IP)== TRUE | IP < 1){
			stop("Data entering error!")
		}
}

I <- IP

###


########
# ============ Define empty vectors (initialisation) ============

max.NB <- array(NA, c(K, N))
NB.t.intrial.std <- NB.t.intrial.int <- NB.t.outtrial.sample	<-  array(NA, c(K, life))
NB.pop.intrial.std <- NB.pop.intrial.int <- NB.pop.outtrial.sample <- rep(NA, K)
NB.pop.sample <- array(NA, c(K, outer.N))

# Initialising variables in order for compiling process to work better
n <- n1 <- n2 <- e1 <- e2 <- mu.D <- tau.D <- mu.HR <- sigma.HR <- size <- j <- 1
HR_PU.X <- rep(NA, inner.N)
inner.X <- matrix(NA, nrow=inner.N, ncol=4) 

# Record of ENPVSI for all lambdas for each sample size.
# Note: first column shows value of lambda for each row. Column names shows number
# of patients in each arm of the trial size (trial size = 2 * that number).
pop.ENPVSI.HR <- array(c(50000, 100000, rep(NA, N * K)), c(K, (N + 1)))
dimnames(pop.ENPVSI.HR)[[2]] <- c("lambda", "200", "400")

# Calculate discounted population NB with current information over life 

# Discounted population NB with current information separately for each year of the lifetime of 
# the technology.
# Note: the choice of treatment is based on current information but in this approach (K.W.'s) 
# for consistency the NB assigned is the one of the posterior NB (as for the outwith not-in-trial 
# patients group), for recording see B9.2. 
# For more information on this approach see comments at B9.2, A2.2 and A3.
NB.pop.eY.current <- array(NA, c(K, life))

# Discounted population NB with current information summed up over whole lifetime of technology.
# Note: = the sum of all 'NB.pop.eY.current' for each outer loop.
NB.pop.current <- array(NA, c(K, outer.N))

# Expected discounted population NB with current information over whole lifetime of technology.
# Note: = the means over all columns of 'NB.pop.current' for each sample size.
ENB.pop.current <- array(NA, c(K, 1, N))

# Expected discounted population NB with sample information over whole lifetime of technology.
# Note: = the means over all columns of 'NB.pop.sample' for each sample size.
ENB.pop.sample <- array(NA, c(K, 1, N))

# The variables below are mostly for checking model functionality

# Record of which treatment is optimal at each lambda, trial size and outer loop.
optimal <- array(c(50000, 100000, rep(NA, N * K)), c(K, N + 1, outer.N))

# Record of mean updated per-patient net benefits for each treatment, lambda, 
# trial size and outer loop.
ENB.IandII <- array(c(50000, 100000, rep(NA, N * K)), c(K, N + 1, 2, outer.N))

# Record of mean in-trial-per-patient net benefits for each treatment, lambda, 
# trial size and outer loop.
ENB.X.IandII <- array(c(50000, 100000, rep(NA, N * K)), c(K, N + 1, 2, outer.N))

# Record of population net benefit of each sample size, outer loop and lambda. 
NB.pop.sample.rec <- array(NA, c(outer.N, N, k))

# Record 'mu.D' for each outer loop of each trial size.
mu.D.rec <- array(NA, c(outer.N, N))

# Record 'mu.D' for each outer loop of each trial size.
# Note: currently not used. 
HR_PU.Act.rec <- array(NA, c(outer.N, N))

# Record real 'trial.time' from the function for each outer loop.
trial.time.rec <- array(NA, c(outer.N, N))


# ============ Compile functions for faster running ============

Run.model.intervention.arm <- cmpfun(Run.model.intervention.arm.uncompiled)
Compute.net.benefit <- cmpfun(Compute.net.benefit.uncompiled)
Run.model.intervention.arm.X <- cmpfun(Run.model.intervention.arm.X.uncompiled)
Compute.net.benefit.X <- cmpfun(Compute.net.benefit.X.uncompiled)
Estimate.trial.time <- cmpfun(Estimate.trial.time.uncompiled)


# ============ Start ENPVSI calculation for parameter HR_PU ============

# Save original parameter draw
saved.HR_PU  <- HR_PU

########
#  B1. Draw sample from prior logHR (='parameter of interest')
# (modified B1.a part 1 in Ades et al. 2004 (page 217) as per P.H.'s model code)
# Note: 1) The hazard ratio is the main parameter of interest in this analysis. 
#		2) The log-normal distribution is applied to a variable whose logarithm
# 		   is normally-distributed. It is usually parameterised with mean and variance, 
#		   or in Bayesian inference, with mean and precision, where precision is the inverse of
#		   the variance.
		tau.0 <- 1 / (HR_PU_SD^2)  # Define precision parameter as: 1/variance of logHR_PU.
		log_HR_PU <- rnorm(outer.N, mu.0, HR_PU_SD)

# (B1.a part 2 in Ades et al. 2004)
# Theta 1 (θ1, here named 't1') drawing:

# Draw a sample baseline parameter theta 1 from its prior distribution.
# Ades et al. 2004 write to use a gamma distribution for hazard rates (t1 is the baseline 
# i.e. standard care event rate) with the expected number of events and total number of years at
# risk as alpha and beta parameters. 
# This would not make sense in our case.
# Here it is instead set to represent the overall (whole model horizon) event rate of 25.1% 
# (expected value for Standard therapy arm) in standard care using a Beta distribution (RSD = 5%).
		H <- Get.BHP.RSD(0.251, 5)
		t1 <- rbeta(outer.N, H[1], H[2])

# (B1.b in Ades et al. 2004)
# Calculate implied prior for theta 2 (here named 't2')
		t2 <- t1 * exp(log_HR_PU)
########


####
# Create progress bar and display that progress bar is for simulations part only 
cat("\nProgress bar for the simulation part of the model:\n")
pb <- txtProgressBar(min = 0, max = N, style = 3)
####


###
# Generate a 'random' number to choose a different random number seed for different instances and computers
systemtime.sec <- format(Sys.time(), "%S")
systemtime.sec <- as.numeric(systemtime.sec)
systemtime.min <- format(Sys.time(), "%M")
systemtime.min <- as.numeric(systemtime.min)
systemtime.hou <- format(Sys.time(), "%H")
systemtime.hou <- as.numeric(systemtime.hou)
seed.number <- systemtime.sec * systemtime.min * systemtime.hou
set.seed(seed.number)
###


# ============ Sample size loop starts here ============
for (n in 1:N) {
#for (n in 1) {

# Define numbers in trial arms for this sample size loop
	n1 <- N1[n]
	n2 <- N2[n]





########
# A1.


# ============ Open outer loop ============
	for (j in 1:outer.N) {  

		########
		# B2. draw sufficient statistic D
		# Note: here this represents the number of patients who get a PU in each arm of the trial.
		e1 <- rpois(1, t1[j] * n1)
		e2 <- rpois(1, t2[j] * n2)
		
		# mu.D = mean of log(HR_PU)
		# Note: 1) The term '(e2 * n1) / (e1 * n2)' is the risk ratio and not the hazard ratio.
		# 		2) Here this represents the risk ratio found in that trial, n1 and n2 being equal cancel each other out. 
		mu.D <- log((e2 * n1) / (e1 * n2))
		
		# Record 'mu.D' for each outer loop of each trial size.
		mu.D.rec[j, n] <- mu.D
		
		# tau.D = precision of log(HR_PU) (increases with number of events e1 and e2)
		tau.D <- 1 / (1 / e1 + 1 / e2)
		########
		
		########
		# B3. and B4.
		
		# B3.1. update prior with sample(D) to give posterior distribution (for parameter of 
		# interest HR_PU)
		HR_PU <- exp(rnorm(inner.N, (mu.0 * tau.0 + mu.D * tau.D) / (tau.0 + tau.D), 
					 sqrt(1/(tau.0 + tau.D))))
		# Note: the SD (i.e. 'sqrt(1 / (tau.0 + tau.D))') decreases with tau.0 and tau.D)
		
		HR_PU.X <- exp(rep(mu.D, inner.N))   # B4.1. mean of DATA
		
		
		# Inner Monte Carlo simulation loop for
		# B3.2 and below for
		# B4.2
# ============ open inner loop ============
		for (i in 1:inner.N) {  
			inner[i, (3:4)] <- Run.model.intervention.arm(i)
			inner.X[i, (3:4)] <- Run.model.intervention.arm.X(i)
		}  
# ============ Close inner loop ============


		# Nothing changes for the results of the standard treatment arm 
		# (compared to base case run from '1_master.R'), therefore they
		# are used here as well.  
		inner.X[, 1] <- inner[, 1]
		inner.X[, 2] <- inner[, 2] 
		
		
		# Calculated net benefits for results of inner Monte Carlo simulation loop for
		# B3.3 and below for
		# B4.3
		NB <- Compute.net.benefit()
		NB.X <- Compute.net.benefit.X()
		########
			
		########
		# B5. record NB of optimum treatment for this outer loop j 
		# (based on updated prior distribution results from B3.)
		# Note: for easier understanding of pop.EVSI results, especially when outer.N is set to 1, 
		# it is included a recording of the optimal strategy in each trial size, the expected NB 
		# using the updated distribution for both strategies (ENB.IandII) and of the expected NB 
		# using the in-trial distribution for both strategies (ENB.X.IandII).
		for (k in 1:K) {
			max.NB[k, n] <- max(mean(NB[, 1, k]), mean(NB[, 2, k]))
			optimal[k, n + 1, j] <- if (mean(NB[, 1, k]) >= mean(NB[, 2, k])){
								1
							}else{
								2
							}
			ENB.IandII[k, n + 1, 1, j] <- mean(NB[, 1, k])
			ENB.IandII[k, n + 1, 2, j] <- mean(NB[, 2, k])
			ENB.X.IandII[k, n + 1, 1, j] <- mean(NB.X[, 1, k])
			ENB.X.IandII[k, n + 1, 2, j] <- mean(NB.X[, 2, k])
		}
		########
		
		
# ============ Population ENPVSI calculations ============		
		
		########
		# B6. Calculate time to trial reporting for this outer loop
		# (rounded up since e.g. the code below would interpret 4.6 as 4 because it only uses 
		# integers and by t=4 the trial would actually not be over yet)
		# Note: After the trial time 12 months could be added in between the trial and 
		# post-reporting phase in the loop below because the function estimates the time-to-event
		# of the trial and not the full 
		# time to trial reporting. It should not be added here since the trial is already over at 
		# that point and all patients are treated with 'current information' optimal alternative.
		
		# Define variables for trial time function and run function
		mu.HR <- mu.D
		sigma.HR <- sqrt(1 / tau.D)
		size <- (N1[n] + N2[n])
		trial.time <- Estimate.trial.time(mu.HR, sigma.HR, size)
		trial.time.rec[j, n] <- trial.time
		trial.time <- ceiling(trial.time / 12)
				
		# Truncate trial time to time of useful lifetime of the intervention, in this cases the
		# trial is not over by the end of the useful lifetime of the intervention technology.
		trial.time <- if (trial.time <= life) trial.time else life
		########
		
		
		########
		# B7, B8 and B9.
		
		# Calculate population NB for this outer loop for the following patient groups:
		# Note: The variable m1 and m2 are used for calculation efficiency and are re-defined/used.
		for (k in 1:K)  {
			
			# per-year calculations
			
			# Patients in trial treated with standard care (intervention 1; group 1) = B7.1
			m1 <- mean(NB.X[, 1, k])
			for (t in 1:life) {
				NB.t.intrial.std[k, t] <- (if (t <= trial.time) {
											m1 * (n1 / trial.time)
										  }else{ 
											0
										  }
										 ) * DF[t]		
			}

			# Patients in trial treated with the intervention (intervention 2; group 2) = B7.2
			m2 <- mean(NB.X[, 2, k])
			for (t in 1:life) {
				NB.t.intrial.int[k, t] <- (if (t <= trial.time) {
											m2 * (n2 / trial.time)
										  }else{
										  0
										  }
										  ) * DF[t]
			}

			# Patients outwith trial (patients not in trial; group 3 and 4) = B8. and B9.1
			m2 <- mean(NB[, 2, k])
			m1 <- mean(NB[, 1, k])
			for (t in 1:life) {
				NB.t.outtrial.sample[k, t] <-   (if (t <= trial.time) {
												# Important note: code here based on premise that 'm1' is for standard therapy.
												    (if (approach == OWR){
														# When using Only-With-Research approach
														if (NBs.current[k, 2] >= NBs.current[k, 1]){
															m2
														}else{
															m1
														}
													 }else{
														# When using Only-In-Research approach
														m1
													 }
													) * (I - n1 / trial.time - n2 / trial.time)
													
												}else{
													max.NB[k, n] * I
											    }
											    ) * DF[t]
			}
			
			# All patients in comparator without trial conducting = B9.2 (step included by K.W.):
			# Record posterior NB for treatment deemed optimal with current (before trials) 
			# information, to use later as comparator for EVSI calculation, see A2.2 and A3.
			# This step was included after finding accordance with Ades et al. 2004: page 222, 
			# section "Computing Strategy" especially last formula and last paragraph.
			#
			# Important note: The treatment deemed optimal with current information CAN be the
			#				  one with the lower NB if OIR approach is chosen (!). 
			# 				  Because an OIR research approach means that currently one has 
			#				  chosen standard therapy for all patients except for those in one 
			#				  arm of the trial. The choice of OIR could be made based on 
			#				  uncertainty on effectiveness of the new treatment or other 
			#				  considerations.
			# Note: code here based on premise that 'm1' is for standard therapy.
			if (approach == OWR){
				# When using Only-With-Research approach
				if (NBs.current[k, 2] >= NBs.current[k, 1]) {
					for (t in 1:life) {
						NB.pop.eY.current[k, t] <- m2 * I * DF[t]
					}
				}else{
					for (t in 1:life) {
						NB.pop.eY.current[k, t] <- m1 * I * DF[t]
					}
				}			
			}else{
				# When using Only-In-Research approach
				for (t in 1:life) {
					NB.pop.eY.current[k, t] <- m1 * I * DF[t]
				}
			}
			# Summing up per-year calculations over lifetime of technology
			
			NB.pop.intrial.std[k] <- sum(NB.t.intrial.std[k, ])
			NB.pop.intrial.int[k] <- sum(NB.t.intrial.int[k, ])
			NB.pop.outtrial.sample[k] <- sum(NB.t.outtrial.sample[k, ])
			
			NB.pop.current[k, j] <- sum(NB.pop.eY.current[k, ])
		}
		########
		
		
		########
		# B10. Record the sum of the expected net benefits over all groups in B7, B8 and B9.
		# Slightly changed from original approach by P.H. in order to be more similar to 
		# the approach of Ades et al. 2004.
		NB.pop.sample[, j] <-  NB.pop.intrial.std + NB.pop.intrial.int + NB.pop.outtrial.sample
		########
		
		# Record of population net benefit of each sample size, outer loop and lambda. 
#		NB.pop.sample.rec[j, n, ] <- NB.pop.sample[, j]
		
	}
# ============ Close outer loop ============
########
			
	
	########
	# A2.
	
	for (k in 1:K) {
	# Both calculations below are set to remove NAs before taking the mean. This was added after
	# one simulation with several outer loops returned NAs as result for one sample size (n = 1477). 
	# This new approach avoids that one NA result in any of the outer loops causes the whole sample
	# size result to be NA.
	
	# A2.1 record population expected NB of a decision based on sample information for this 'n' 	
		ENB.pop.sample[k, 1, n] <- mean(NB.pop.sample[k, ], na.rm = TRUE)
	# A2.2 record population expected NB of a decision based on current information
	# Only necessary with K.W.'s approach. It uses strategies based on current information but 
	# NB of this strategy based on (used from) mean NB with updated information - for 
	# consistency & comparability.
		ENB.pop.current[k, 1, n] <- mean(NB.pop.current[k, ], na.rm = TRUE)
	}
	########

	########
	# A3. Record population ENPVSI for sample size loop n
	# As reference use the treatment alternative from current information with updated 
	# NB values for this alternative in order to compare equal entities at the same times.
	pop.ENPVSI.HR[, n + 1] <- ENB.pop.sample[, , n] - ENB.pop.current[, , n]
	########
	
	
	# Update progress bar
    setTxtProgressBar(pb, n)
	
	# Save backup of results in a file for this outer loop.
	write.csv(pop.ENPVSI.HR, 
			  file=paste(as.character(n), ".csv", sep=""), row.names=FALSE)
}
# ============ Sample size loop ends here ============


# Save and display 'pop.ENPVSI.HR'
write.csv(pop.ENPVSI.HR, file=paste("pop_ENPVSI_HR.csv", sep=""), row.names=FALSE)
cat("\n\nThe object open outer loop was saved \nin the file 'pop_ENPVSI_HR.csv'.\n")
cat("\nThe values in 'pop.ENPVSI.HR' are:\n")
print(pop.ENPVSI.HR)



# Restore original parameter draw 
HR_PU  <- saved.HR_PU

# Save some other results in separate files
source('save_other_e_results.R')


# ============ Final timing part ============  

# Tell the R software to display only 2 decimals.
options(digits=2) 

# Calculate running time of script by subtracting processor time at the start
# from current processor time.
duration <- matrix(c(start.time - proc.time()))

# Display the time it took to run the script 
cat("\n \n Running this script took ", c( - duration[3, ]), "seconds ")
cat("(~", c(( - duration[3, ]) / 60), "minutes). \n")

options(digits=5)  # Reset decimals to 5

