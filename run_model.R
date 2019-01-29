
# Name of script: run_model.R
# Purpose:	To run a discrete-state and time Markov model for two 
#			alternative strategies that produces discounted cost & benefit 
#			outcomes (but without drawing of parameter values or evaluation of 
#			outputs/outcomes - for that see other scripts).
# Author(s):	Peter Hall, Klemens Wallner
# Note(s):	1) The script is split into a (small) constant part that needs to 
#			run only once before the first simulation and a (large) function 
#			part that needs to run every time the 'run model' function is 
#			called.
#			2) Some error detection sections have been added. This includes
#			the checking for illogical values in the transition probabilities.


# ============ Constant parts of the script ============
# Things the function does not need to do every time it is called.
# (Definitions that are completely (!) independent of the current 
# simulation run (value of 'i')).


# List of states with index numbers:

# State No. 1 ..... Hospital PU-free
# State No. 2 ..... Hospital with one PU
# State No. 3 ..... Hospital with two PUs
# State No. 4 ..... Hospital PU-free, after having had at least one PU
# State No. 5 ..... Discharged with one PU
# State No. 6 ..... Discharged with two PUs
# State No. 7 ..... At home PU-free
# State No. 8 ..... Dead


########
# Transition matrices

# Reference list for transition probabilities (here from state 1)
# tps[1, 1, t] <- 0
# tps[1, 2, t] <- 0
# tps[1, 3, t] <- 0
# tps[1, 4, t] <- 0
# tps[1, 5, t] <- 0
# tps[1, 6, t] <- 0
# tps[1, 7, t] <- 0
# tps[1, 8, t] <- 0

# The probability to stay in one state (tps[S, S, t]) is listed at the very 
# end of the calculation loop for easy calculation and should not be 
# listed before. 

# Initialise variables 
# Note that the transition probabilities are set to 0 so that later
# in each simulation only the ones that are not 0 have to be defined.
tps <- array(0, c(S, S, T))
tpi <- array(0, c(S, S, T))
mr.t <- rep(NA, T)
traces <- matrix(0, nrow = T, ncol = S)
tracei <- matrix(0, nrow = T, ncol = S)

# Define mortality depending on starting age and model year (cycle)
for (t in 1:T) {
	mr.t[t] <- 	mr[startage + ceiling(t/cyclesperyear)]
}

# Define Dead state as absorbing state
tps[8, 8, 1:T] <- 1
tpi[8, 8, 1:T] <- 1
########


########
# Define discount factor array and discount factor matrix.
DF <- c(1/((1 + disc)^(0:199)))  # Time horizon, up to 200 years.
# Discount factor matrix for 'cyclesperyear' cycles per year and a 
# time horizon of 'T/cyclesperyear' years
dm <- diag(c(rep(DF[1:(T/cyclesperyear)], each=cyclesperyear)))
########


########
# Calculate the discounted life expectancy (LE) at age = start age + 3 years

# Create discount factor matrix 
# Note: The discounting up to t = 3 years is already done in the 
# 		discounted traces.
DF.matrix <- matrix(c(DF), nrow=1, ncol=length(DF))
#(1 X 200)

# Create array were each year of the remaining LE can 
# be multiplied with a different discount factor
# E.g. for 3.41 years = (1,1,1,1,0.41)
LE_sep <- c(rep(1,floor(LE.at.markov.end)), (LE.at.markov.end - floor(LE.at.markov.end)))

# Variable for making code line below it more clear
# Note: It calculates the number or elements (i.e. length) of the 
#		DF vector that are not used for discounting the LE (i.e. are 0).
#		This is needed in order to flexibly create an object that
# 		can be matrix multiplied with a DF row matrix.
zeros <- length(DF) - length(LE_sep)

# Create LE matrix that matches the DF matrix
LE.matrix_match <- matrix(c(LE_sep, rep(0, zeros)), nrow=length(DF), ncol=1)
#(200 X 1)

# Calculate discounted LE
dLE <- (DF.matrix %*% LE.matrix_match)
########


# ============ Function part of the script ============
# What to do every time the function is called. 
# Note: the function is uncompiled and will be compiled, i.e. translated
# into hardware language, before running it)


Run.model.uncompiled <- function(i) {

########
# Transition matrix - standard care arm

# Cycle-independent transition probabilities

tps[1, 2, 1:T] <- (1 - pDis[i]) * pPU[i] * (1 - pPU2[i]) * 2
tps[1, 3, 1:T] <- (1 - pDis[i]) * pPU[i] * pPU2[i]
tps[1, 5, 1:T] <- pDis[i] * pPU[i] * (1 - pPU2[i]) * 2
tps[1, 6, 1:T] <- pDis[i] * pPU[i] * pPU2[i]
tps[1, 7, 1:T] <- pDis[i] * (1 - pPU[i]) * (1 - pPU2[i])

tps[2, 3, 1:T] <- (1 - pDis[i]) * (1 - pHeal[i]) * pPU2[i]
tps[2, 4, 1:T] <- (1 - pDis[i]) * pHeal[i] * (1 - pPU2[i])
tps[2, 5, 1:T] <- pDis[i] * (pHeal[i] * pPU2[i] + (1 - pHeal[i]) * (1 - pPU2[i]))
tps[2, 6, 1:T] <- pDis[i] * (1 - pHeal[i]) * pPU2[i]
tps[2, 7, 1:T] <- pDis[i] * pHeal[i] * (1 - pPU2[i])

tps[3, 2, 1:T] <- (1 - pDis[i]) * (1 - pHeal[i]) * pHeal[i] * 2
tps[3, 4, 1:T] <- (1 - pDis[i]) * pHeal[i] * pHeal[i]
tps[3, 5, 1:T] <- pDis[i] * (1 - pHeal[i]) * pHeal[i] * 2
tps[3, 6, 1:T] <- pDis[i] * (1 - pHeal[i]) * (1 - pHeal[i])
tps[3, 7, 1:T] <- pDis[i] * pHeal[i] * pHeal[i]

tps[4, 2, 1:T] <- (1 - pDis[i]) * pRec[i] * (1 - pPU2[i]) * 2
tps[4, 3, 1:T] <- (1 - pDis[i]) * pRec[i] * pPU2[i]
tps[4, 5, 1:T] <- pDis[i] * pRec[i] * (1 - pPU2[i]) * 2
tps[4, 6, 1:T] <- pDis[i] * pRec[i] * pPU2[i]
tps[4, 7, 1:T] <- pDis[i] * (1 - pRec[i]) * (1 - pPU2[i])

tps[5, 2, 1:T] <- pRA[i] * (1 - pPU2D[i])
tps[5, 3, 1:T] <- pRA[i] * pPU2D[i]
tps[5, 6, 1:T] <- (1 - pRA[i]) * (1 - min(pHeal[i] * HR_heal[i],1)) * pPU2D[i]
tps[5, 7, 1:T] <- (1 - pRA[i]) * min(pHeal[i] * HR_heal[i],1) * (1 - pPU2D[i])

tps[6, 2, 1:T] <- pRA[i] * min(pHeal[i] * HR_heal[i],1)
tps[6, 3, 1:T] <- pRA[i] * (1 - min(pHeal[i] * HR_heal[i],1))
tps[6, 5, 1:T] <- (1 - pRA[i]) * min(pHeal[i] * HR_heal[i],1) * (1 - min(pHeal[i] * HR_heal[i],1)) * 2
tps[6, 7, 1:T] <- (1 - pRA[i]) * min(pHeal[i] * HR_heal[i],1) * min(pHeal[i] * HR_heal[i],1)


# Cycle-dependent transition probabilities

###
# Define mortality in hospital and discharged states
for (t in 1:time_ICU_ST[i]){
tps[1:4, 8, t] <- mICU[i]
}

aux <- time_ICU_ST[i] + 1
for (t in aux:T){
tps[1:4, 8, t] <- mr.t[t]
}

for (t in 1:T){
tps[5:7, 8, t] <- mr.t[t]
}

tps[c(2, 3, 5, 6), 8, 1:T] <- tps[c(2, 3, 5, 6), 8, 1:T] * HR_D.PU[i]
tps[c(3, 6), 8, 1:T] <- tps[c(3, 6), 8, 1:T] * HR_D.PU[i]
###

# Ensuring that all elements in transition matrices are between 0 and 1. For this we use the 
# parallel minimum and parallel maximum function (which compares element row individually).
tps <- pmin(tps, 1)
tps <- pmax(tps, 0)

for (t in 1:T) {
	for (s in 1:7) {
		tps[s, s, t] <- max(1 - sum(tps[s, -s, t]),0)
	}
}
########

########
# Transition matrix - intervention arm

# Cycle-independent transition probabilities

tpi[1, 2, 1:T] <- (1 - pDis[i]) * min(pPU[i] * HR_PU[i],1) * (1 - min(pPU2[i] * HR_PU[i],1)) * 2
tpi[1, 3, 1:T] <- (1 - pDis[i]) * min(pPU[i] * HR_PU[i],1) * min(pPU2[i] * HR_PU[i],1)
tpi[1, 5, 1:T] <- pDis[i] * min(pPU[i] * HR_PU[i],1) * (1 - min(pPU2[i] * HR_PU[i],1)) * 2
tpi[1, 6, 1:T] <- pDis[i] * min(pPU[i] * HR_PU[i],1) * min(pPU2[i] * HR_PU[i],1)
tpi[1, 7, 1:T] <- pDis[i] * (1 - min(pPU[i] * HR_PU[i],1)) * (1 - min(pPU2[i] * HR_PU[i],1))

tpi[2, 3, 1:T] <- (1 - pDis[i]) * (1 - pHeal[i]) * min(pPU2[i] * HR_PU[i],1)
tpi[2, 4, 1:T] <- (1 - pDis[i]) * pHeal[i] * (1 - min(pPU2[i] * HR_PU[i],1))
tpi[2, 5, 1:T] <- pDis[i] * (pHeal[i] * min(pPU2[i] * HR_PU[i],1) + (1 - pHeal[i]) * (1 - min(pPU2[i] * HR_PU[i],1)))
tpi[2, 6, 1:T] <- pDis[i] * (1 - pHeal[i]) * pPU2[i]
tpi[2, 7, 1:T] <- pDis[i] * pHeal[i] * (1 - pPU2[i])

tpi[3, 2, 1:T] <- (1 - pDis[i]) * (1 - pHeal[i]) * pHeal[i] * 2
tpi[3, 4, 1:T] <- (1 - pDis[i]) * pHeal[i] * pHeal[i]
tpi[3, 5, 1:T] <- pDis[i] * (1 - pHeal[i]) * pHeal[i] * 2
tpi[3, 6, 1:T] <- pDis[i] * (1 - pHeal[i]) * (1 - pHeal[i])
tpi[3, 7, 1:T] <- pDis[i] * pHeal[i] * pHeal[i]

tpi[4, 2, 1:T] <- (1 - pDis[i]) * pRec[i] * (1 - pPU2[i]) * 2
tpi[4, 3, 1:T] <- (1 - pDis[i]) * pRec[i] * pPU2[i]
tpi[4, 5, 1:T] <- pDis[i] * pRec[i] * (1 - pPU2[i]) * 2
tpi[4, 6, 1:T] <- pDis[i] * pRec[i] * pPU2[i]
tpi[4, 7, 1:T] <- pDis[i] * (1 - pRec[i]) * (1 - pPU2[i])

tpi[5, 2, 1:T] <- pRA[i] * (1 - pPU2D[i])
tpi[5, 3, 1:T] <- pRA[i] * pPU2D[i]
tpi[5, 6, 1:T] <- (1 - pRA[i]) * (1 - min(pHeal[i] * HR_heal[i],1)) * pPU2D[i]
tpi[5, 7, 1:T] <- (1 - pRA[i]) * min(pHeal[i] * HR_heal[i],1) * (1 - pPU2D[i])

tpi[6, 2, 1:T] <- pRA[i] * min(pHeal[i] * HR_heal[i],1)
tpi[6, 3, 1:T] <- pRA[i] * (1 - min(pHeal[i] * HR_heal[i],1))
tpi[6, 5, 1:T] <- (1 - pRA[i]) * min(pHeal[i] * HR_heal[i],1) * (1 - min(pHeal[i] * HR_heal[i],1)) * 2
tpi[6, 7, 1:T] <- (1 - pRA[i]) * min(pHeal[i] * HR_heal[i],1) * min(pHeal[i] * HR_heal[i],1)


# Cycle-dependent transition probabilities

###
# Define mortality in hospital and discharged states
for (t in 1:time_ICU_IN[i]){
tpi[1:4, 8, t] <- mICU[i] * HR_D.SeP[i]
}

aux <- time_ICU_IN[i] + 1
for (t in aux:T){
tpi[1:4, 8, t] <- mr.t[t] * HR_D.SeP[i]
}


for (t in 1:T){
tpi[5:7, 8, t] <- mr.t[t]
}

tpi[c(2, 3, 5, 6), 8, 1:T] <- tpi[c(2, 3, 5, 6), 8, 1:T] * HR_D.PU[i]
tpi[c(3, 6), 8, 1:T] <- tpi[c(3, 6), 8, 1:T] * HR_D.PU[i]
###

# Ensuring that all elements in transition matrices are between 0 and 1. For this we use the 
# parallel minimum and parallel maximum function (which compares element row individually).
tpi <- pmin(tpi, 1)
tpi <- pmax(tpi, 0)


for (t in 1:T) {
	for (s in 1:7) {
		tpi[s, s, t] <- max(1 - sum(tpi[s, -s, t]),0)
	}
}
########


########
# Check transition probabilities for illogical values.

####
# Check that all transition matrix rows add up to 1, if this is not the 
# case collect and display information on the error and force the rows 
# to sum up to 1.

# Create vectors containing the row sums over all rows of
# the TP array. 
# Note: 1) Below R calculates the sums of each row over the whole 
#		model horizon (3rd dimension). Therefore each of the 
#		sums should add up to the number of cycles of the model horizon i.e. 
#		'T' (if they add up to 1 within each cycle).
# 		2) The number sums = number of states = S.
# 		3) Calculations in R are usually correct even beyond 10 digits after 
#		the dot. However, to avoid artificial inaccuracies, which we found in
#		test runs we are here rounding to 12 digits.

	aux_s <- round(rowSums(tps),12)
	aux_i <- round(rowSums(tpi),12)

	if ((length(aux_s[aux_s != T]) > 0) || (length(aux_i[aux_i != T]) > 0)){ 
		rowSumsError <- 1

# Re-scaling mechanism - force all TP rows to add up to 1.
		for (t in 1:T) {
			for (s in 1:S) {
				tps[s, , t] <- tps[s, , t] / sum(tps[s, , t])
				tpi[s, , t] <- tpi[s, , t] / sum(tpi[s, , t])
			}
		}
	}else{
		rowSumsError <- 0
	}
####

########


# ============ CREATE MARKOV TRACES ============

# Markov trace standard care ('traces') and intervention ('tracei') arm

# Set trace for t = 1
# Set starting proportion for each non-empty state 
# i.e. distribution of starting population as percentages (1% = 0.01)
# (The starting population for all other states remains at 0.)
traces[1, 1] <- tracei[1, 1] <- 1 - rPU
traces[1, 2] <- tracei[1, 2] <- rPU


# Calculate trace for t >= 2
# Note: this code block could possibly be made faster through vectorisation, however
# 		so far I did not find a way to vectorise it. 
for (t in 2:T) {
	TM.S <- tps[, , (t - 1)]
	traces[t, ] <- traces[t - 1, ] %*% TM.S
	TM.I <- tpi[, , (t - 1)]
	tracei[t, ] <- tracei[t - 1, ] %*% TM.I
}


# ============ Outputs ============

########
# QALYs up to t = 3y

# For models with constant, utilities for each state and each transition reward.

# Create discounted trace (T X S) for each treatment alternative
dtraces <- dm %*% traces
dtracei <- dm %*% tracei



# Get state utilities in form of column matrices
# Note: Currently the utilities are the same for both treatment strategies. 
# Alternatively one can define the row of the object to be used: e.g. uS1[1, i] and uS1[2, i]).
U.s <- matrix(c(uS1[i], uS2[i], uS3[i], uS4[i], uS5[i], uS6[i], uS7[i], uS8[i]), 
			  nrow=S, ncol=1)
U.i <- matrix(c(uS1[i], uS2[i], uS3[i], uS4[i], uS5[i], uS6[i], uS7[i], uS8[i]), 
			  nrow=S, ncol=1)

# Transforming annual utilities into per-cycle utilities.
U.s <- U.s / cyclesperyear
U.i <- U.i / cyclesperyear

# Note:
# For getting the QALYs for each state and each cycle use a diagonal 
# matrix (S x S) with the utility vector as the diagonal.
#QALY.c.s <- dtraces %*% diag(U.s)

###
# Half-cycle corrected traces
# Note: Currently only the first cycle is corrected because 
# the last one is not in the markov model. 

# Discounted traces
HCC_dtraces <- dtraces
HCC_dtraces[1,] <- 0.5 * HCC_dtraces[1,]
HCC_dtracei <- dtracei
HCC_dtracei[1,] <- 0.5 * HCC_dtracei[1,]

# Undiscounted traces
HCC_traces <- dtraces
HCC_traces[1,] <- 0.5 * HCC_traces[1,]
HCC_tracei <- dtracei
HCC_tracei[1,] <- 0.5 * HCC_tracei[1,]
###

# Calculate QALYs from states 
QALY.s <- sum(HCC_dtraces %*% U.s)
QALY.i <- sum(HCC_dtracei %*% U.i)

# Adding transition utilities
# Note: 1) For transition utilities that only occur from specific states the 
#		calculation is a bit more difficult. E.g. only for transitions from 
#		state 1 to state 2:
#		QALY.s <- QALY.s + (uS1S2[i] * (matrix(dtraces[1:(T - 1), 1], nrow=1)
#					 %*% matrix(tps[1, 2, 1:(T - 1)], ncol=1)))
QALY.s <- QALY.s + (sum(diff(dtraces[, 8])) * utoS8[i])
QALY.i <- QALY.i + (sum(diff(dtracei[, 8])) * utoS8[i])
########



########
# Life years up to t = 3y
# un-discounted life years
LY.s <- (sum(HCC_traces[, -8]) / cyclesperyear)
LY.i <- (sum(HCC_tracei[, -8]) / cyclesperyear)

#IncLY <- LY.i - LY.s

# discounted life years
#dLY.s <- (sum(HCC_dtraces[, - 8]) * cyclesperyear)
#dLY.i <- (sum(HCC_dtracei[, - 8]) * cyclesperyear)
########


########
# Costs up to t = 3y

# Get state costs in form of column matrices
# Note: states 1 to 4 use ICU costs for all cycles, correction see below. 
C.s <- matrix(c(cS1[1, i], cS2[1, i], cS3[1, i], cS4[1, i], cS5[1, i], cS6[1, i], cS7[1, i], cS8[1, i]), 
			  nrow=S, ncol=1)
C.i <- matrix(c(cS1[2, i], cS2[2, i], cS3[2, i], cS4[2, i], cS5[2, i], cS6[2, i], cS7[2, i], cS8[2, i]), 
			  nrow=S, ncol=1)
  
# Calculate costs for all states (assuming cycle independence)
# Note: 1) Here the cycle independence means (lower) non-ICU costs used throughout. This is corrected below.
#       2) This is a work around solution for easier coding.
COST.s <- sum(HCC_dtraces %*% C.s)
COST.i <- sum(HCC_dtracei %*% C.i)


###
# Adjust costs for ICU cycles having higher costs
# For those cycles only add the difference of the drawn costs for SHB and ICU bed. 

# Create traces that are empty beyond the cycle in ICU 

HCC_dtraces.aux <- HCC_dtraces
HCC_dtraces.aux[((time_ICU_ST[i]+1):T),] <- 0 

HCC_dtracei.aux <- HCC_dtracei
HCC_dtracei.aux[((time_ICU_IN[i]+1):T),] <- 0 

# Calculate additional per cycle costs from patients being in ICU
c_diff <- cICU[i] - cSHB[i]

# Create cost vector for costs difference only (same for each strategy)
C.both <- matrix(0, nrow=S, ncol=1)
C.both[1:4] <- c_diff

# Calculate additional costs for whole traces
COST.s.add <- sum(HCC_dtraces.aux %*% C.both)
COST.i.add <- sum(HCC_dtracei.aux %*% C.both)

# Add up cycle independent and additional costs to real total costs
COST.s <- COST.s + COST.s.add
COST.i <- COST.i + COST.i.add
###

# Additional comments:
# For getting the costs for each state and each cycle separately use a diagonal
# matrix (S x S) with the cost vector as the diagonal.
#Costs.c.s <- dtraces %*% diag(C.s)


# Adding transition cost (currently not used!)
# Note: 1) For transition cost that only occur from specific states the calculation is a bit 
#		more difficult. e.g. only for transitions from state 1 to state 2:
#		COST.s <- COST.s + (cS1S2[i] * (matrix(dtraces[1:(T - 1), 1], nrow=1) %*% 
#							matrix(tps[1, 2, 1:(T - 1)], ncol=1)))

#COST.s <- COST.s + (sum(diff(dtraces[, 8])) * ctoS8[i])
#COST.i <- COST.i + (sum(diff(dtracei[, 8])) * ctoS8[i])
########


########
# Life-time expansion - add QALYs and costs for remaining life time for 
# each arm (after t = 3y) according to distribution in states.
# We take population distribution and assume they live for their 
# remaining LE in the same states.

###
# QALYs after t = 3y

# Transforming per-cycle utilities into annual utilities.
U.s_a <- U.s * cyclesperyear
U.i_a <- U.i * cyclesperyear

# Record of discounted utilities (annual value)
dU.at3y.s <- dtraces[T,] %*% U.s_a
dU.at3y.i <- dtracei[T,] %*% U.i_a

# Calculate remaining QALYs after the original model horizon
QALY.s.rem <- dU.at3y.s * dLE
QALY.i.rem <- dU.at3y.i * dLE

# Calculate QALYs for total model horizon 
QALY.total.i <- QALY.i + QALY.i.rem
QALY.total.s <- QALY.s + QALY.s.rem
###

###
# Life years after t = 3y
# un-discounted life years

traces[T,] * LE.at.markov.end
tracei[T,] * LE.at.markov.end

LY.s <- LY.s + (sum(traces[T, -8]) * LE.at.markov.end)
LY.i <- LY.i + (sum(tracei[T, -8]) * LE.at.markov.end)
###

###
# Costs after t = 3y

# Transforming per-cycle costs into annual costs.
C.s_a <- C.s * cyclesperyear
C.i_a <- C.i * cyclesperyear

# Record of discounted costs (annual value)
dC.at3y.s <- dtraces[T,] %*% C.s_a
dC.at3y.i <- dtracei[T,] %*% C.i_a

# Calculate remaining costs after the original model horizon
COST.s.rem <- dC.at3y.s * dLE
COST.i.rem <- dC.at3y.i * dLE

# Calculate COSTs for total model horizon 
COST.total.s <- COST.s + COST.s.rem
COST.total.i <- COST.i + COST.i.rem
###

########


########
# Calculate model calibration and check variables

cv <- rep(NA, 8)

if (evppi.run != 1){
# Calculation mean duration in Hospital (in days, 1, 2, 3, 4; Standard Care).
cv[1] <- cyclelength * sum(traces[, c(1, 2, 3, 4)])

# Calculation mean duration of being with PU (in days, 2, 3, 5, 6; any # of PUs; standard care arm).
# Note: This is per patient value over all patients (sic!) although 
# not all those patients have had a PU. 
cv[2] <- cyclelength * sum(traces[, c(2, 3, 5, 6)])

# Currently not used.
cv[3] <- 0


# Calculation of % of patients that do develop at least one new PU when 
# having none at start (standard care)
cv[4] <- (sum(tps[1, c(2, 3, 5, 6), T])) * (sum(traces[, 1]))

# Calculation of % of patients that do develop at least one new PU when 
# having none at start (intervention)
cv[5] <- (sum(tpi[1, c(2, 3, 5, 6), T])) * (sum(tracei[, 1]))


# Calculation of % of patients that do develop a 2nd PU (from anywhere) (standard)
cv[6] <- sum((sum(tps[1, c(3, 6), T])) * (sum(traces[, 1])), 
				   (sum(tps[2, c(3, 6), T])) * (sum(traces[, 2])), 
				   (sum(tps[4, c(3, 6), T])) * (sum(traces[, 4])), 
				   (sum(tps[5, c(3, 6), T])) * (sum(traces[, 5])))

# Calculation of % of patients that do develop a 2nd PU (only from PU2 or A0) (standard)
cv[7] <- sum((sum(tps[2, c(3, 6), T])) * (sum(traces[, 2])), 
		  (sum(tps[4, c(3, 6), T])) * (sum(traces[, 4])), (sum(tps[5, c(3, 6), T])) * (sum(traces[, 5])))

# Error dummy - set to 1 (0) for an error (no error) in this simulation
cv[8] <- rowSumsError
}
########

model_return <- c(QALY.total.s, COST.total.s, QALY.total.i, COST.total.i, LY.s, LY.i, cv, COST.s.rem, COST.i.rem, QALY.s.rem, QALY.i.rem)

# Model return:
return(model_return)
}


########
# The part below is an unused part form 
# Peter Hall's script 'model.R' and might be useful to keep.
# It is however commented out now and not updated with the current 
# parameter names. 

#IncCOST <- COST.i - COST.s

#### 
# ICER
#ICER.LYs <- IncCOST/IncLYs
#ICER.QALYs <- IncCOST/IncQALYs
####


# Diagnostic returns:

#return(array(c(LYs.s, LYs.i, IncLYs, QALYs.s, QALYs.i, IncQALYs, COST.s, COST.i, 
#				IncCOST, ICER.LYs, ICER.QALYs, "NA"), c(3, 4)))

#return(tracei[is.na(tracei)==TRUE], )
#return(TM.I[is.na(TM.I)==TRUE])

#return(c(QALYs.s, COST.s, QALYs.i, COST.i, C.i, U.i))
########


