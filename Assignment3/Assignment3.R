setwd("D:/Masters Program Chalmers/Projects and Labs/CMBS")

# Q1 - d ------------------------------------------------------------------

dayCount <- c(162,267,271,185,111,61,27,8,3,1)
deathNotices <- c(0:9)

functionG <- function(p, lambda1, lambda2){
  num <- p*(lambda1^deathNotices)*exp(-lambda1)
  denom <- num + (1-p)*(lambda2^deathNotices)*exp(-lambda2)
  functionG <- num/denom
}

N <- 10000
pOld <- 0.5
lambda1Old <- rbeta(1,1,1)
lambda2Old <- rbeta(1,1,1)

PVector <- rep(0, times = N)
lambda1Vector <- rep(0, times = N)
lambda2Vector <- rep(0, times = N)

for (iter in 1:N){
  numerator <- sum(dayCount*functionG(pOld,lambda1Old, lambda2Old))
  p <- numerator/sum(dayCount)
  lambda1 <- sum(deathNotices*dayCount*functionG(pOld,lambda1Old, lambda2Old))
  lambda1 <- lambda1/(1 + numerator)
  lambda2 <- sum(deathNotices * dayCount * (1 - functionG(pOld,lambda1Old, lambda2Old)))
  lambda2 <- lambda2/(1 + sum(dayCount*(1 - functionG(pOld,lambda1Old, lambda2Old))))
  
  pOld <- p
  PVector[iter] <- p
  lambda1Old <- lambda1
  lambda1Vector[iter] <- lambda1
  lambda2Old <- lambda2
  lambda2Vector[iter] <- lambda2
}

#Prediction using above p, lambda_1 and lambda_2
predictedCounts <- sum(dayCount)*(pOld*dpois(deathNotices, lambda1Old) + 
  (1 - pOld)*dpois(deathNotices, lambda2Old))
legend(x = 5.5, y= 225, legend = c("Observed", "Predicted by bipoisson"),col = c("black","blue"), pch = c(20,20))

#Plotting predicted counts, observed counts and new predicted counts if there was one distribution
# i.e. Assignment 2 Q2-b

# If there was only one distribution
lambda <- sum(deathNotices*dayCount)/sum(dayCount)
newPredictedCounts <- sum(dayCount)*dpois(deathNotices,lambda)
plot(deathNotices, y = dayCount, col = "black", pch = 20, ylim = c(0,300),
     cex = 1.5,xlab = "Death notices", ylab = "Frequency", main = "Observed and Expected counts")
lines(deathNotices, y = predictedCounts, col = "blue", pch = 20, cex = 1.5)
lines(deathNotices, y = newPredictedCounts, col = "red", pch = 20, cex = 1.5)
legend(x = 5.5, y= 225, legend = c("Observed","Predicted by bipoisson","Predicted by unipoisson"),col = c("black","blue","red"), pch = c(20), lty = c(0,1,1))


differenceCount <- dayCount - predictedCounts
differenceCount

par(mfrow = c(2,2))
plot(PVector, ylab = "p", main = "p Value convergence")
plot(lambda1Vector, ylab = expression(paste(lambda, "1")), main = expression(paste(lambda, "1 convergence")))
plot(lambda2Vector, ylab = expression(paste(lambda, "2")), main = expression(paste(lambda, "2 convergence")))
plot(deathNotices,differenceCount,type = "b", ylab = "Difference", xlab = "Number of deaths",main = "Observed - Predicted")
dev.off()
