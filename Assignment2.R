setwd("D:/Masters Program Chalmers/Projects and Labs/CMBS")

# Assignment 2 ------------------------------------------------------------

# Q1 - a ------------------------------------------------------------------

alphaParameter = 1.5
betaParameter = 0.1
quantiles = c(13, 30)
probQ1a = pgamma(quantiles[2], alphaParameter, betaParameter) - pgamma(quantiles[1], alphaParameter, betaParameter)
probQ1a

# Q1 - b ------------------------------------------------------------------

frequencyVelocity <- c(10,100,1000,10000,100000)
windVelocities = vector(mode = "list", length = length(frequencyVelocity))
windPower = vector(mode = "list", length = length(frequencyVelocity))
expectedWindPower <- rep(0, length(frequencyVelocity))

windPowerFunction <- function(velocity){
  if(velocity >= 13 && velocity < 23){
    wP <- cos(velocity*pi/10 - 3*pi/10) + 1
    return(wP)
  }
  else if(velocity >= 23 && velocity < 30){
    wP <- (1031 + 46*velocity - velocity^2)/780
    return(wP)
  }
  else{
    return(0)
  }
}
set.seed(777)
for (iFreq in 1:length(frequencyVelocity)){
  windVelocities[[iFreq]] <- rgamma(frequencyVelocity[iFreq], alphaParameter, betaParameter)
  windPower[[iFreq]] <- rep(0,times = frequencyVelocity[[iFreq]])
  windPower[[iFreq]] <- sapply(windVelocities[[iFreq]], windPowerFunction)
  expectedWindPower[iFreq] <- mean(windPower[[iFreq]])
}
expectedWindPower
#print expected wind power

# Bootstrap for CI
nBootSamples <- 1000
bootStrapSamples <- vector(mode = "list", length = length(frequencyVelocity))
bootStrapSampleMean <- vector(mode = "list", length = length(frequencyVelocity))
diffMean <- vector(mode = "list", length = length(frequencyVelocity))
quantiles <- vector(mode = "list", length = length(frequencyVelocity))
confidenceInterval <- vector(mode = "list", length = length(frequencyVelocity))
set.seed(777)
for (iFreq in 1:length(frequencyVelocity)) {
  sampleValues <- sample(windPower[[iFreq]], nBootSamples*frequencyVelocity[iFreq], replace = TRUE)
  bootStrapSamples[[iFreq]] <- matrix(sampleValues, nrow = frequencyVelocity[iFreq], ncol = nBootSamples)
  bootStrapSampleMean[[iFreq]] <- colMeans(bootStrapSamples[[iFreq]])
  diffMean[[iFreq]] <- expectedWindPower[iFreq] - bootStrapSampleMean[[iFreq]]
  quantiles[[iFreq]] <- quantile(diffMean[[iFreq]],c(0.25,0.975))
  confidenceInterval[[iFreq]] <- expectedWindPower[iFreq] - c(quantiles[[iFreq]][2], quantiles[[iFreq]][1])
}
#print confidence interval


# Q1 - c ------------------------------------------------------------------

velocitySequence <- seq(0,60,length.out = 1000)
densityVelocity <- dgamma(velocitySequence,alphaParameter, betaParameter)
plot(velocitySequence,densityVelocity, ylab = "Density", xlab = "velocity", main = "Prior density")
powerSequence <- sapply(velocitySequence, windPowerFunction)
plot(velocitySequence,powerSequence, ylab = "Wind power", xlab = "velocity", main = "Likelihood")
posteriorSequence <- powerSequence*densityVelocity
plot(velocitySequence, posteriorSequence,, ylab = "Density", xlab = "velocity", main = "Posterior upto some constant")
areaUnderTheCurve <- sum(posteriorSequence) * 60 / 1000 
areaUnderTheCurve

# Q1 - d ------------------------------------------------------------------


windPowerIntegrateFn <- function(velocity){
  if(velocity >= 13 && velocity < 23){
    wP <- (cos(velocity*pi/10 - 3*pi/10) + 1) * (betaParameter^alphaParameter)*velocity^(alphaParameter - 1)*exp(-betaParameter*velocity)/gamma(alphaParameter)
    return(wP)
  }
  else if(velocity >= 23 && velocity < 30){
    wP <- ((1031 + 46*velocity - velocity^2)/780) * (betaParameter^alphaParameter)*velocity^(alphaParameter - 1)*exp(-betaParameter*velocity)/gamma(alphaParameter)
    return(wP)
  }
  else{
    return(0)
  }
} 
integrate(windPowerIntegrateFn,lower = 13, upper = 30)


# Q1 - e ------------------------------------------------------------------

windPowerPosterior <- sapply(velocitySequence, windPowerIntegrateFn)
plot(velocitySequence, windPowerPosterior, ylab = "Expected wind power upto a constant", xlab = "velocity", main = "Posterior and proposal fucntion")
points(velocitySequence, dgamma(velocitySequence, 469/25, 23/25)/3, col = "red")
legend(50, 0.03, legend=c("Posterior", "Proprosal"),
       col=c("black", "red"), pch = 1, cex=0.8)

# Q1 - f ------------------------------------------------------------------

alphaImportance <- 469/25
betaImportance <- 23/25
importanceSamplingFn <- function(velocity){
  if(velocity >= 13 && velocity < 23){
    wP <- (cos(velocity*pi/10 - 3*pi/10) + 1) * (betaParameter^alphaParameter)*velocity^(alphaParameter - 1)*exp(-betaParameter*velocity)/gamma(alphaParameter)
    wP <- wP/((betaImportance^alphaImportance)*velocity^(alphaImportance - 1)*exp(-betaImportance*velocity)/gamma(alphaImportance))
    return(wP)
  }
  else if(velocity >= 23 && velocity < 30){
    wP <- ((1031 + 46*velocity - velocity^2)/780) * (betaParameter^alphaParameter)*velocity^(alphaParameter - 1)*exp(-betaParameter*velocity)/gamma(alphaParameter)
    wP <- wP/((betaImportance^alphaImportance)*velocity^(alphaImportance - 1)*exp(-betaImportance*velocity)/gamma(alphaImportance))
    return(wP)
  }
  else{
    return(0)
  }
} 
sampleImportance <- rgamma(10000,alphaImportance, betaImportance)
imporanceFunctionValue <- sapply(sampleImportance, importanceSamplingFn)
mean(imporanceFunctionValue)


# Q2 - b ------------------------------------------------------------------

dayCount <- c(162,267,271,185,111,61,27,8,3,1)
deathNotices <- c(0:9)
lambda <- sum(deathNotices*dayCount)/sum(dayCount)
predictedCounts <- sum(dayCount)*dpois(deathNotices,lambda)
predictedCounts
plot(deathNotices, y = dayCount, col = "black", pch = 20,cex = 1.5,xlab = "Death notices", ylab = "Frequency", main = "Observed and Expected counts")
points(deathNotices, y = predictedCounts, col = "red", pch = 20, cex = 1.5)
legend(x = 7.5, y= 225, legend = c("Observed", "Predicted"),col = c("black","red"), pch = c(20,20))


# Q2 - g ------------------------------------------------------------------

lengthPrior <- 1
N <- 100

posteriorLambda1 <- rep(0, times = N)
posteriorLambda2 <- rep(0, times = N)
p <- rep(0, times = N)

posteriorLambda1[1] <- rgamma(lengthPrior,1,1)
posteriorLambda2[1] <- rgamma(lengthPrior,1,1)
p[1] <- runif(lengthPrior)

Zs <- function(Y,i, p, l1, l2){
  numerator =   p*(l1^i)*exp(-l1)
  denom = p*(l1^i)*exp(-l1) + (1-p)*(l2^i)*exp(-l2)
  prob = numerator/denom
  return(rbinom(1, Y, prob))
}

estimatedZ <- rep(0,times = 10)
iVector <- c(0:9)

for (i in 2:N){
  for (iZ in 1:length(deathNotices)){
    estimatedZ[i] <- Zs(dayCount[iZ], iZ - 1, p[i-1],posteriorLambda1[i - 1], posteriorLambda2[i - 1] )
  }
  p[i] <- rbeta(1, sum(estimatedZ) + 1, sum(dayCount - estimatedZ) + 1)
  posteriorLambda1[i] <- rgamma(1,sum(iVector*estimatedZ) + 1, sum(estimatedZ) + 1)
  posteriorLambda2[i] <- rgamma(1,sum(iVector*estimatedZ) + 1, sum(estimatedZ) + 1)
}
hist(p)
hist(posteriorLambda1, ylab =  "Frequency", xlab = "Parameter value", main = "Histogram for lambda_1")
hist(posteriorLambda2, ylab =  "Frequency", xlab = "Parameter value", main = "Histogram for lambda_2")

plot(posteriorLambda1[1:10000], ylab = "Parameter value", xlab = "Iteration number", main = "Trace plot for first 10000 iterations for lambda_1")
plot(posteriorLambda2[1:10000], ylab = "Parameter value", xlab = "Iteration number", main = "Trace plot for first 10000 iterations for lambda_2")
plot(p[1:100], ylab = "Parameter value", xlab = "Iteration number", main = "Trace plot for first 100 iterations for p")

avgP <- mean(p[2001:N])
lambda1 <- mean(posteriorLambda1[2001:N])
lambda2 <- mean(posteriorLambda2[2001:N])
avgP
lambda1
lambda2
