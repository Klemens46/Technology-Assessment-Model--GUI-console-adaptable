
# Name of script: run_model.intervention.arm.R
# Purpose:	To run the intervention part of a discrete-state and time Markov 
#			model for two alternative strategies that produces discounted 
#			cost & benefit outcomes (but without drawing of parameter values
#			or evaluation of outputs/outcomes - for that see other scripts).
# Author(s):	Peter Hall, Klemens Wallner
# Note(s):	1) The script is split into a (small) constant part that needs to 
#			run only once before the first simulation and a (large) function 
#			part that needs to run every time the 'run model' function is 
#			called.
#			2) Some error detection sections have been added. This includes
#			the checking for illogical values in the transition probabilities.


# List of states with index numbers
# State No. 1 ..... Hospital PU-free or with stage 1 PU(s)
# State No. 2 ..... Hospital with stage 2 or above PU
# State No. 3 ..... Hospital with two stage 2 or above PU
# State No. 4 ..... Hospital PU-free or with stage 1 PU(s), after having had at least one PU
# State No. 5 ..... Discharged with stage 2 or above PU
# State No. 6 ..... Discharged with two stage 2 or above PU
# State No. 7 ..... At home PU-free
# State No. 8 ..... Dead


# ============ TRANSITION MATRICES ============

# Reference list for transition probabilities (here from state 1)
# tps[1, 1, t] <- 0
# tps[1, 2, t] <- 0
# tps[1, 3, t] <- 0
# tps[1, 4, t] <- 0
# tps[1, 5, t] <- 0
# tps[1, 6, t] <- 0
# tps[1, 7, t] <- 0
# tps[1, 8, t] <- 0

# The probability to stay in one state (tps[S, S, t]) is listed at the very end of 
# the calculation loop for easy calculation and should not also be listed before (duplication). 



# ============ Function part of the script ============
# What to do every time the function is called. 
# Note: the function is uncompiled and will be compiled, i.e. translated
# into hardware language, before running it)

Run.model.intervention.arm.X.uncompiled <- function(i) {

########
# Transition matrix - intervention arm

# Cycle-independent transition probabilities

tpi[1, 2, 1:T] <- (1 - pDis[i]) * min(pPU[i] * HR_PU.X[i],1) * (1 - min(pPU2[i] * HR_PU.X[i],1)) * 2
tpi[1, 3, 1:T] <- (1 - pDis[i]) * min(pPU[i] * HR_PU.X[i],1) * min(pPU2[i] * HR_PU.X[i],1)
tpi[1, 5, 1:T] <- pDis[i] * min(pPU[i] * HR_PU.X[i],1) * (1 - min(pPU2[i] * HR_PU.X[i],1)) * 2
tpi[1, 6, 1:T] <- pDis[i] * min(pPU[i] * HR_PU.X[i],1) * min(pPU2[i] * HR_PU.X[i],1)
tpi[1, 7, 1:T] <- pDis[i] * (1 - min(pPU[i] * HR_PU.X[i],1)) * (1 - min(pPU2[i] * HR_PU.X[i],1))

tpi[2, 3, 1:T] <- (1 - pDis[i]) * (1 - pHeal[i]) * min(pPU2[i] * HR_PU.X[i],1)
tpi[2, 4, 1:T] <- (1 - pDis[i]) * pHeal[i] * (1 - min(pPU2[i] * HR_PU.X[i],1))
tpi[2, 5, 1:T] <- pDis[i] * (pHeal[i] * min(pPU2[i] * HR_PU.X[i],1) + (1 - pHeal[i]) * (1 - min(pPU2[i] * HR_PU.X[i],1)))
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

for (t in 1:T){
	tpi[1:7, 8, t] <- mr.t[t]
}

tpi[c(2, 3, 5, 6), 8, 1:T] <- tps[c(2, 3, 5, 6), 8, 1:T] * HR_D.PU[i]
tpi[c(3, 6), 8, 1:T] <- tps[c(3, 6), 8, 1:T] * HR_D.PU[i]


# Ensuring that all elements in the transition matrices are between 0 and 1. For this we use the 
# parallel minimum and parallel maximum function (which compares element row individually).
pmin(tpi, 1)
pmax(tpi, 0)


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

	aux_i <- round(rowSums(tpi),12)

	if (length(aux_i[aux_i != T]) > 0){ 
		rowSumsError <- 1

# Re-scaling mechanism - force all TP rows to add up to 1.
		for (t in 1:T) {
			for (s in 1:S) {
				tpi[s, , t] <- tpi[s, , t] / sum(tpi[s, , t])
			}
		}
	}else{
		rowSumsError <- 0
	}
####

########


# ============ MARKOV TRACES ============

# Markov trace intervention ('tracei') arm

# Set trace for t = 1
# Set starting proportion for each non-empty state 
# i.e. distribution of starting population as percentages (1% = 0.01)
# (The starting population for all other states remains at 0.)
tracei[1, 1] <- 1 - rPU
tracei[1, 2] <- rPU


# Calculate trace for t >= 2
for (t in 2:T) {
	TM.I <- tpi[, , (t - 1)]
	tracei[t, ] <- tracei[t - 1, ] %*% TM.I
}


# ============ Outputs ============

########
# QALYs up to t = 3y

# For models with constant, utilities for each state and each transition reward.

# Create discounted trace (T X S) for each treatment alternative
dtracei <- dm %*% tracei

# Get state utilities in form of column matrices
# Note: Currently the utilities are the same for both treatment strategies. 
# Alternatively one can define the row of the object to be used: e.g. uS1[1, i] and uS1[2, i]).
U.i <- matrix(c(uS1[i], uS2[i], uS3[i], uS4[i], uS5[i], uS6[i], uS7[i], uS8[i]), 
			  nrow=S, ncol=1)

# Transforming annual utilities into per-cycle utilities.
U.i <- U.i / cyclesperyear

# Calculate utilities from states 
QALY.i <- sum(dtracei %*% U.i)


########
# Costs up to t = 3y

# Get state costs in form of column matrices
C.i <- matrix(c(cS1[2, i], cS2[i], cS3[i], cS4[i], cS5[i], cS6[i], cS7[i], cS8[i]), 
			  nrow=S, ncol=1)

# Calculate cost from being in any of the states
COST.i <- sum(dtracei %*% C.i)
########


########
# Life-time expansion - add QALYs and costs for remaining life time for 
# each arm (after t = 3y) according to distribution in states.
# We take population distribution and assume they live for their 
# remaining LE in the same states.

###
# QALYs after t = 3y

# Transforming per-cycle utilities into annual utilities.
U.i_a <- U.i * cyclesperyear

# Record of discounted utilities (annual value)
dU.at3y.i <- dtracei[T,] %*% U.i_a

# Calculate remaining QALYs after the original model horizon
QALY.i.rem <- dU.at3y.i * dLE

# Calculate QALYs for total model horizon 
QALY.i <- QALY.i + QALY.i.rem
###


###
# Costs after t = 3y

# Transforming per-cycle costs into annual costs.
C.i_a <- C.i * cyclesperyear

# Record of discounted costs (annual value)
dC.at3y.i <- dtracei[T,] %*% C.i_a

# Calculate remaining costs after the original model horizon
COST.i.rem <- dC.at3y.i * dLE

# Calculate COSTs for total model horizon 
COST.i <- COST.i + COST.i.rem
###

########


# Model return:
return(c(QALY.i, COST.i))
}

