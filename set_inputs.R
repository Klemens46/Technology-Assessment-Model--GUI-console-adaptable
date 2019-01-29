
# Name of script: set_inputs.R
# Purpose:	Defining and drawing input values for the stochastic model parameters from 
#			distributions for use in the model in the script 'run_model.R' (some 
#			non-stochastic values are also defined here).
# Author(s):	Peter Hall, Klemens Wallner
# Note(s):	1) For defining transition probabilities there is also a deterministic
#			version, which in probabilistic runs (if 'prob.run == 1') just provides the 
#			distribution means and is then overwritten.
#			2) For background mortality there is only a deterministic version.


# ============ TRANSITION PROBABILITIES & LIFE EXPECTANCY ============

# Define ratio (proportion) of patients that have a PU at start
rPU <- 0

########
# Defining background mortality
# Note: only 2nd table is used, 1st one is for reference purposes

###
# Setting up life table (annual): non-weighted mean from male and female life table values from age 0 to
# 110 i.e. annual probability to die at a certain age, 
# e.g. mr[1] = annual probability to die before the first birthday
# 		mr[81] = annual probability when being 80 to die before the 81st birthday.
# Note: value reported at age 110 was 1 therefore which might help with model validation.

# Read in mean Canadian life table values (annual)
filepath <- c("C:\\Users\\User3\\Dropbox\\R Software\\PU CEA model\\Alberta PU model\\MR and LE over age for m&f in Canada.csv")
D1 <- read.csv(file=paste(filepath, sep = ""))  # Here resulting in a 'data.frame' type object
mr_a <- D1$mr_a  # Calls the element named 'mr' in the list 'D1'
###


# Setting up life table (daily - using same value for whole year; calculated from 'mr_a' above): 
# non-weighted mean from male and female life table values from age 0 to 
# 130 i.e. daily probability to die at a certain age, 
# e.g. mr[1] = daily probability to die before the first birthday;
#		mr[81] = daily probability when being 80 to die before the 81st birthday.
# Note: values are rounded to speed up calculations
mr <- 1 - ((1 - mr_a)^(1/cyclesperyear))
mr <- round(mr, 10)
########



########
# Life expectancy (LE) in years for male and females at a certain age. E.g. the first value is the
# LE at age 0 (i.e. at birth). The 5th value is the LE at the 4th birthday etc.  

# Read in mean Canadian life expectancy (LE) values in years
# Note: Currently uses same file as above.
LE <- D1$LE  # Calls the element named 'LE' in the list 'D1'
# Note: There are 111 elements in this vector.
########


# Life expectancy (LE; in years) after model horizon i.e. at start age 
# ('startage') plus 3 years (= 'T'). Calculated using 'LE' but with adding 4 because
# the LE variable starts at age 0. E.g. 3 years after being 0 one is 3 
# and the value at 3 is the 4th value in 'LE' (value at 0, 1, 2, 3).
LE.at.markov.end <- LE[startage + 4]
# Note: Works up to a start age of 107 years.


########
# In a deterministic model run define probability parameters (also set means for probabilistic run):

# Probability of being discharged i.e. probability that initial condition is sufficiently treated.
pDis <- 0.182

# Probability of being re-admitted to hospital i.e. that a PU (the PUs) is (are) considered too 
# difficult to treat outside of a hospital.
# Currently not possible i.e. 0.
pRA <- 0

# Probability of developing a first PU in Hospital.
pPU <- 0.037100

# Probability of developing a second PU in Hospital.
pPU2 <- 0.070155
#pPU2 <- 0.07997113

# Probability of developing a second PU when discharged with one PU.
# Currently assumed to be 0 because they are assumed to be mobile 
# again after treatment of their initial condition. 
pPU2D <- 0

# Hazard ratio of developing a PU in intervention arm 'HR_PU'. The pre-trial estimated 
# PU incidence in ST and IT arms (25.1% vs. 16.817%) results in a risk ratio of 0.67 (16.817/25.1).
# Note: A risk ratio is different from a hazard ratio: a hazard ratio is applied each time 
# a hazard is encountered (here each cycle) - a risk ratio compares the cumulative risks of 
# two groups over a period of time (here the full time horizon of the Markov model).
mu.0 <- log(0.592559)
log_HR_PU <- mu.0
HR_PU <- exp(log_HR_PU)

# Probability of PU healing.
pHeal <- 0.034236424

# Hazard ratio of PU healing when discharged with a PU.
log_HR_heal <- log(1)
HR_heal <- exp(log_HR_heal)

# Probability of PU recurrence while being PU free and not discharged as PU free.
# Note: (Currently the same as pPU).
pRec <- pPU * 1

# Time in ICU for patients staying longer than 2 cycles
# Note: 1)This is the time after which patients that are in hospital states are 
# 		 considered to be transferred to non-ICU wards.
#		2) Currently it is set to a minimum of 6 days.
#		3) Currently this is set to be drawn from a Poisson distribution with mean of 6.
#		4) Since we have not data on this we assume equality in both treatment arms.
time_ICU_ST <- time_ICU_IN <- pmax(rpois(inner.N, 6), 6)
time_ICU_ST <- time_ICU_IN <- round(time_ICU_ST/cyclelength)


# Mortality  in ICU without PU.
mICU <- 0.09997 # per 3 day cycle

# Hazard ratio of death due to the IES device (also as proxy for all hypothetical side effects)
log_HR_D.SeP <- log(1)
HR_D.SeP <- exp(log_HR_D.SeP)

# Hazard ratio of death when having a PU.
log_HR_D.PU <- log(1.92)
HR_D.PU <- exp(log_HR_D.PU)
########



# Re-define depending on model setting:
if (prob.run == 1) {
########
# In a probabilistic model run:
#
# Note: 'Get.BHP.RSD' is a custom function that converts a 
# mean and a given RSD into the corresponding Beta distribution 
# hyper-parameters (see script 'get.beta.hyperparameters.R').
# 
# Caution: The 'Get.BHP.RSD' might not work when the mean is 0 or 1.

# Probability of being discharged i.e. probability that initial condition is 
# sufficiently treated and therefore patient is discharged from hospital.
H <- Get.BHP.RSD(pDis, 5)
pDis <- rbeta(inner.N, H[1], H[2])

# Probability of PU healing.
H <- Get.BHP.RSD(pHeal, 10)
pHeal <- rbeta(inner.N, H[1], H[2])

# Probability of developing a first PU in hospital.
H <- Get.BHP.RSD(pPU, 10)
pPU <- rbeta(inner.N, H[1], H[2])

# Hazard ratio of developing a PU in intervention arm 'HR_PU'. 
# The pre-trial estimated PU incidence in ST and IT arms (25.1% vs. 16.817%) results in a risk 
# ratio of 0.67 (16.817/25.1).
# Notes: 1) A risk ratio is different from a hazard ratio: a hazard ratio is applied each time 
# 			a hazard is encountered (here each cycle) - a risk ratio compares the cumulative risks of 
# 			two groups over a period of time (here the full time horizon of the Markov model).
# 		 2)	mean = 0.59; RSD = 50%.
HR_PU_SD <- 0.29
log_HR_PU <- rnorm(inner.N, mu.0, HR_PU_SD)
HR_PU <- exp(log_HR_PU)
#HR_PU <- rlnorm(inner.N, - 0.0.286254, 0.246221)


# Probability of developing a second PU in hospital.
H <- Get.BHP.RSD(pPU2, 10)
pPU2 <- rbeta(inner.N, H[1], H[2])

# Probability of being re-admitted to hospital i.e. that a PU (the PUs) is (are)
# considered too difficult to treat outside of a hospital.
# Currently not possible i.e. 0.
pRA <- rep(0, inner.N)

# Probability of developing a second PU when discharged with one PU.
# Currently assumed to be 0 because they are assumed to be mobile 
# again after treatment of their initial condition. 
pPU2D <- rep(0, inner.N)

# Hazard ratio of PU healing when discharged with a PU.
# (mean = 0.95; RSD = 5%)
log_HR_heal <- rnorm(inner.N, log(0.95), 0.049969)
HR_heal <- exp(log_HR_heal)

# Probability of PU recurrence while being PU free and not discharged as PU free.
# (Currently the same as pPU).
pRec <-  pPU * 1

# Mortality  in ICU without PU.
H <- Get.BHP.RSD(mICU, 10)
mICU <- rbeta(inner.N, H[1], H[2])

###
# Hazard ratio of death due to the IES device
# Note: 1) this is also used as proxy for all hypothetical side effects)
#		2) mean = 1; RSD = 5%
log_HR_D.SeP <- rnorm(inner.N, log(1), 0.05)
HR_D.SeP <- exp(log_HR_D.SeP)
###

# Hazard ratio of death when having a PU.
# (RSD is calculated from 95% CI (non-symmetric) and 'Goal-Seek' tool in Excel).
# (mean = 1.92; RSD = 12.46%)
log_HR_D.PU <- rnorm(inner.N, log(1.92), 0.08128)
HR_D.PU <- exp(log_HR_D.PU)
########
}



# ============ UTILITY PARAMETERS ============

# Notes:	1) Some of the code below is only usable if the utilities for each state stay 
# 			constant within each model run (for all cycles).
# 			2) Notation: uS1 means the utility of the state one; 
#			uS1S2 means a utility as a transition reward when transiting from state 1 to 2; 
#  			utoS8 means the same transition reward for every transition to state 8.
#  			u = utility, H = in hospital, D = discharged from hospital, A = at home; 
#  			PU means one PU, PU2 means two PUs. 
#  			uHa = in hospital PU free after having at least one PU; 
# 			3) The utilities mostly have two variable names with identical
#			values, for easier identification of where the utility applies.
# 			4) The state utilities are the same in the standard care and intervention arms.
# 			5) Highest utility state: At home disease free without PU. From this value
# 			the other disutility increments subtracted.
# 			6) IMPORTANT: drawing needs to be in order of highest utility first
#			and then descending!


# List of states with index numbers (ordered from highest utility to lowest):
# State No. 7 ..... At home PU-free
# State No. 1 ..... Hospital PU-free
# State No. 4 ..... Hospital PU-free, after having had at least one PU
# State No. 5 ..... Discharged with  PU
# State No. 2 ..... Hospital with  PU
# State No. 6 ..... Discharged with two PUs
# State No. 3 ..... Hospital with two PUs
# State No. 8 ..... Dead


########
# Pre-calculations for utilities

####
# Defining distribution means

# Assumption of HRQoL utility at starting age of 52 years: 0.90
u52 <- 0.9

# Assumption disutility from being in hospital: 0.05
uDisutH <- 0.05

# Disutility from having only one PU = 0.2903 (utility = 0.7097)
uDisutPU2 <- 0.2903

# Disutility from having two PUs = 0.43545 (utility = 0.56455): 
# full disutility of first and half of second.
uDisutPU22 <- uDisutPU2 + 0.5 * uDisutPU2
####

# Mean utility from hospital and PU status times baseline age-dependent utility
S1 <- S4 <- u52 * (1 - uDisutH)
S2 <- 		u52 * (1 - uDisutH - uDisutPU2) 
S3 <- 		u52 * (1 - uDisutH - uDisutPU22)  
S7 <- 		u52 * (1)
S5 <- 		u52 * (1 - uDisutPU2)   
S6 <- 		u52 * (1 - uDisutPU22)  

####
# Defining and utility decrements based on utility differences between states
# Note: 'Get.GHP.RSD' is a custom function that converts a 
# mean and a given RSD into the corresponding Gamma distribution 
# hyper-parameters (see script 'get.gamma.hyperparameters.R').

# Defining relative standard deviation (in %) for the decrements
RSDm <- 5
 
uDec.S7.S1.mean <- S7 - S1
H <- Get.GHP.RSD(uDec.S7.S1.mean, RSDm)  
uDec.S7.S1 <- rgamma(inner.N, H[1], H[2])

uDec.S4.S5.mean <- S4 - S5
H <- Get.GHP.RSD(uDec.S4.S5.mean, RSDm)
uDec.S4.S5 <- rgamma(inner.N, H[1], H[2])

uDec.S5.S2.mean <- S5 - S2
H <- Get.GHP.RSD(uDec.S5.S2.mean, RSDm)
uDec.S5.S2 <- rgamma(inner.N, H[1], H[2])

uDec.S2.S6.mean <- S2 - S6
H <- Get.GHP.RSD(uDec.S2.S6.mean, RSDm)
uDec.S2.S6 <- rgamma(inner.N, H[1], H[2])

uDec.S6.S3.mean <- S6 - S3
H <- Get.GHP.RSD(uDec.S6.S3.mean, RSDm)
uDec.S6.S3 <- rgamma(inner.N, H[1], H[2])
####

########

###
# State utilities

# Get hyperparameters for state with highest utility
H <- Get.BHP.RSD(u52, 2.5)

# Calculate for each state
uS7 <- uA0 <- rbeta(inner.N, H[1], H[2])
uS1 <- uH0 <- uS7 - uDec.S7.S1
uS4 <- uHa <- uS1 
uS5 <- uD2 <- uS4 - uDec.S4.S5
uS2 <- uH2 <- uS5 - uDec.S5.S2
uS6 <- uD22 <- uS2 - uDec.S2.S6
uS3 <- uH22 <- uS6 - uDec.S6.S3
uS8 <- rep(0, inner.N)
###


# Transition utilities

# Disutility typical at the end of life
uTerminal.mean <- 0.2  # positive value of it
H <- Get.GHP.RSD(uTerminal.mean, RSDm)
utoS8 <- - rgamma(inner.N, H[1], H[2]) # negative value of it



# ============ COST PARAMETERS ============

########
# Bed
# Daily hospital bed costs:

# For ICU bed per day
cPD_ICU <- 4224.63
 
# For standard hospital bed
cPD_SHB <- 2189



# Equivalent costs per cycle:
cPC_ICU <- cPD_ICU * cyclelength
cPC_SHB <- cPD_SHB * cyclelength



# RSD for standard hospital bed costs (in percent)
RSDm <- 20

# Probabilistic version of the per-cycle costs above
H <- Get.GHP.RSD(cPC_ICU, RSDm)
cICU <- rgamma(inner.N, H[1], H[2])
H <- Get.GHP.RSD(cPC_SHB, RSDm)
cSHB <- rgamma(inner.N, H[1], H[2])
########


########
# Mattress
# Daily costs for IES device (additional to normal ICU costs):

# Mean cost per day (depending on durability assumptions)
cPD_INT <- cPD_SeP <- 33  # Based on daily cost value supplied by trial team

# Mean cost per day (depending on durability assumptions)
cPD_STT <- 0  # Currently this cost variable is not used but kept for future use. 
#cPD_STT <- 0.050432733 # Based on a durability of 9 years


# Equivalent costs per cycle:
cPC_INT <- cPD_INT * cyclelength

# RSD for mattress costs (in percent)
RSDm <- 15

# Probabilistic version of the per-cycle costs above
cSTT <- rep(0, inner.N)
H <- Get.GHP.RSD(cPC_INT, RSDm)
cINT <- rgamma(inner.N, H[1], H[2])
########


########
# PU-specific treatment costs
# Weighted excess mean treatment costs per patient in 
# 2014 BPD for patients with either one or two PU(s) of the described average severity.

####
# In hospital

# Per day with one PU
cHPD <- 10041 / 13.94
# Note: Weighted average of costs from Pham et al. 2011 with adjustment divided by average model LOS in hospital

# Equivalent cost per cycle:
cHPC <- cHPD * cyclelength

# Probabilistic version of the per-cycle costs above
H <- Get.GHP.RSD(cHPC, 15)
# (one PU = 'cH2' and two PUs = 'cH22')
cH <- rgamma(inner.N, H[1], H[2])

# Per day with two PUs
# Currently assumed to be 150% of cost with one PU
cH2 <- cH * 1.5
####





####
# Once discharged from hospital
# Note: presumed equal to hospital costs due to transport etc.
cD <- cH
cD2 <- cH2
####
########


########
# State costs

# Pre-calculations
cS1_STT <- cSHB + cSTT
cS1_INT <- cSHB + cINT

cS2_STT <- cS1_STT + cH
cS2_INT <- cS1_INT + cH

cS3_STT <- cS1_STT + cH2
cS3_INT <- cS1_INT + cH2

cS5_STT <- cD
cS5_INT <- cD

cS6_STT <- cD2
cS6_INT <- cD2

# Actual state costs
# The state costs might be different in the standard care and intervention arms. 
# In those cases the first row represents the standard care and the second the intervention costs.
cS1 <- cH0 <- matrix(c(cS1_STT, cS1_INT), nrow=2, ncol=inner.N, byrow=TRUE)
cS2 <- cH2 <- matrix(c(cS2_STT, cS2_INT), nrow=2, ncol=inner.N, byrow=TRUE)
cS3 <- cH22 <- matrix(c(cS3_STT, cS3_INT), nrow=2, ncol=inner.N, byrow=TRUE)
cS4 <- cHa <- cS1
cS5 <- cD2 <-  matrix(c(cS5_STT, cS5_INT), nrow=2, ncol=inner.N, byrow=TRUE)
cS6 <- cD22 <- matrix(c(cS6_STT, cS6_INT), nrow=2, ncol=inner.N, byrow=TRUE)
cS7 <- cA0 <- matrix(0, nrow=2, ncol=inner.N, byrow=TRUE)
cS8 <- cS7
########


########
# Transition costs (currently not used!)

# Extra costs typical at the end of life
# Note: Converted and inflation adjusted from PRESSURE2 trial model
#cTerminal.mean <- 33614.81
#H <- Get.GHP.RSD(cTerminal.mean, 5)
#ctoS8 <- cTerminal <- rgamma(inner.N, H[1], H[2])
########

