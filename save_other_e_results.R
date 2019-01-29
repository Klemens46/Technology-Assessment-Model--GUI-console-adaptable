
# Name of script: save_other_e_results.R 
# [to run type or copy: source('save_other_e_results.R')]
# Purpose: Save several result objects from running the script '2_compute_enpvsi.R' as CSV files.
# Author(s):	Klemens Wallner


min1oneoptimal <- rep(NA, outer.N)

for (j in 1:outer.N){
min1oneoptimal[j] <- 	if((sum(optimal[, 2, j])==4) && (sum(optimal[, 3, j])==4)){
			0
		}else{
			1
		}
}

# Calculate % of trials where optimal choice was not the expected one from base case,
# for each lambda (K) and each sample size (N).
notasENBbasecase <- array(NA, c(2,2))
not <- array(c(20000, 50000, rep(NA, N * K)), c(K, N + 1, outer.N))

for (k in 1:K){
	aux <- NBs.current[k, 3]

	for (n in 1:N){
		for (j in 1:outer.N){
			not[k, n+1, j] <- 	if((optimal[k, n+1, j])==aux){
							0
						}else{
							1
						}
		}
	notasENBbasecase[k,n] <- mean(not[k,n+1,])
	}
}


x <- c(mean(uS7), mean(uS1), mean(uS4), mean(uS5), mean(uS2), mean(uS6), mean(uS3), mean(uS8))

cat("\n\nThe objects 'min1oneoptimal' and 'x' (mean utility values ranked) were created.\n\n")

options(digits=5)
cat("\n\nThe mean of 'min1oneoptimal' is ", c(mean(min1oneoptimal)), ".\n\n")

cat("\n\nThe % of trials where optimal choice was not the expected one from base case is:\n")

print(notasENBbasecase)
cat("\n\n")


write.csv(optimal, file=paste("optimal.csv", sep=""), row.names=FALSE)
write.csv(trial.time.rec, file=paste("trial.time.rec.csv", sep=""), row.names=FALSE)
write.csv(max.NB, file=paste("max.NB.csv", sep=""), row.names=FALSE)
write.csv(mu.D.rec, file=paste("mu.D.rec.csv", sep=""), row.names=FALSE)
write.csv(log_HR_PU, file=paste("log_HR_PU.csv", sep=""), row.names=FALSE)
write.csv(t1, file=paste("t1.csv", sep=""), row.names=FALSE)
write.csv(t2, file=paste("t2.csv", sep=""), row.names=FALSE)
write.csv(min1oneoptimal, file=paste("min1oneoptimal.csv", sep=""), row.names=FALSE)
write.csv(notasENBbasecase, file=paste("notasENBbasecase.csv", sep=""), row.names=FALSE)

cat("\n\nThe objects 'optimal', 'trial.time.rec', 'max.NB', 'mu.D.rec', \n'log_HR_PU', 't1', 't2', 'min1oneoptimal' and 'notasENBbasecase' were saved as CSV files.\n\n")



