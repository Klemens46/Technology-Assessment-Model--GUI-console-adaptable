
# Name of script: estimate_trial_time.R
# Purpose:	Stochastic forecasting of time-to-event in trial (Markov model based) to be used in 
# 			the script '2_compute_enpvsi.R'.
# Author(s):	Peter Hall, Klemens Wallner


# ============ Constant parts of the script ============

########
# Model overview

# State transition model for stochastic forecasting of time-to-event trial

# 6 states:
# A = patients to be recruited
# B = patients in control arm
# C = patients in intervention arm
# D = patients lost to follow-up
# E = patients having experienced an event
# F = patients that have left the hospital without either getting a PU or being lost to follow-up 
#     (are considered no longer at risk of event) 

# Transition Matrix
#		A					B			C			D		E		F
# A 	1 - sum(r[1:(t - 1)])	r[t]/2		r[t]/2		0		0		0
# B 	0					0			0			f		p		1 - f - p
# C 	0					0			0			f		p.HR	1 - f - p * HR
# D 	0					0			0			1		0		0
# E 	0					0			0			0		1		0
# F		0					0			0			0		0		1

########


########
#
#
# Define model structure
TH <-  40			# time horizon (months)
II <-  100			# number of simulations

## Define input parameters
# Calculated to get p to have a mean of 0.251 (updated baseline STT probability of getting a 
# PU over whole model horizon)
alp <- 251
bet <- 749

# Simplification: If patients get a new PU they get it within one month of entering the 
# trial (ICU patients rarely stay longer than one month in ICU).

e <- c(21, 42) # required number of events 
# Note:  In absence of this in the current trial protocol, it was calculated as expected mean incidence times trial size.
#e <- c(300, 445) # required number of events (Source: Table 1 in protocol, for trial sizes of 902 and 1996 respectively)
f <- 0.005			# monthly probability of being lost to follow-up

# Create empty objects
#h <- NA			# h = monthly hazard of event
p <- NA				# p = monthly probability of event
HR <- NA			# HR = hazard ratio

# Markov trace (separately for each state)
A <- rep(0, TH)
B <- rep(0, TH)
C <- rep(0, TH)
D <- rep(0, TH)
E <- rep(0, TH)
F <- rep(0, TH)
E.abs <- rep(NA, TH)
tau <- rep(NA, II)
########


# ============ Trial time function start  ============
Estimate.trial.time.uncompiled <- function(mu.HR, sigma.HR, size) {

# Assumed monthly recruitment
# Note: Currently calculated as expected sample i.e. 200 divided by expected trial duration (i.e. 15 months)
ar <- 13


# r = monthly recruitment vector
# (Number of patients (in multiples of the assumed monthly recruitment 'ar') that need to be 
#  recruited in every month to add up to the whole trial size. Each month either the assumed 
#  monthly recruitment number of patients is recruited if more patients are needed or 0 when they are not.)
r <- c(rep(ar, ceiling(size / ar)), rep(0, TH - (floor(size / ar))))


# Markov trace (A, B, C, D, E, F)

# A = patients to be recruited
# (At t = 1 nobody has yet been recruited so the whole number of patients in 'r' still need to be recruited.)
A[1] <- sum(r)

## Run simulation 
for (i in 1:II) {

	p <- rbeta(1, alp, bet)
	HR <- rlnorm(1, mu.HR, sigma.HR) 

	for (t in 2:TH) {
# 6 states:
# A = patients to be recruited
# B = patients in control arm
# C = patients in intervention arm
# D = patients lost to follow-up
# E = patients having experienced an event
# F = patients that have left the hospital without either getting a PU or being lost to follow-up 
#     (are considered no longer at risk of event) 
		A[t] <- A[1] - sum(r[1:(t - 1)])
		B[t] <- r[t] / 2
		C[t] <- r[t] / 2
		D[t] <- D[t - 1] + (B[t - 1] + C[t - 1]) * f
		E[t] <- E[t - 1] + (B[t - 1] + C[t - 1] * HR) * p
		F[t] <- F[t - 1] + B[t - 1] * (1 - f - p) + C[t - 1] * (1 - f - p * HR)
}

	E <- round(E, 0)
		
# Create a vector with number of additional patients with event that are still needed in any month.
	EN <- (e[n] - E)

# Create a vector with the time in months it takes to have the desired amount of events.
# Note: When the number of events is not reached within the maximum trial period then the result is 151.
	tau[i] <- length(EN[EN > 0]) + 1 
}

return(mean(tau))
}

# Currently unused code:
#plot(density(tau, bw=1), lty=1, lwd=2, main = "Sampled distribution for Tau", xlab="time (months)")

