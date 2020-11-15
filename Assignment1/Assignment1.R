
# Q2a ---------------------------------------------------------------------

theta <- seq(1.51,1.53,length = 100)
mu <- c(1.5163, 1.5197, 1.5203)
sigmaA <- c(0.001,0.001,0.005)
sigmaL <- 0.001
weights <- c(0.33,0.57,0.1)
densityValues <- function(theta, mu, sigma, weights){
  dValues = rep(0,length(theta))
  for (i in 1:3){
    dValues <- dValues + weights[i]*dnorm(theta, mu[i], sigma[i])
    dValues
  }
  return(dValues)
}
plot(theta,densityValues(theta,mu,sigmaA,weights), xlab = "RI", ylab = expression(pi(theta)), main = "Prior distribution")


# Q2b ---------------------------------------------------------------------

sigmaB <- c(0.000002,0.000002,0.000026)
sigmaB <- sqrt(sigmaB)
plot(theta,densityValues(theta,mu,sigmaB,weights), xlab = "RI", ylab = expression(pi(x)),main = "Prior predictive distribution")


# Q2c ---------------------------------------------------------------------

denominator = 0
x = 1.52083
for (i in 1:3){
  A = (weights[i])/sqrt((sigmaA[i])^2 + sigmaL^2)
  B = exp((-0.5*((x-mu[i])^2))/((sigmaA[i])^2 + sigmaL^2))
  C = A*B
  denominator = denominator + C
}
{numerator = c(0,0,0)
D = c(0,0,0)
E = D
G = E
H = G
weightsQc = numerator}
for (i in 1:3){
  D[i] = x^2 * sigmaL^2 + mu[i]^2*sigmaA[i]^2
  E[i] = (x * sigmaL^2 + mu[i]*sigmaA[i]^2)^2
  E[i] = E[i]/(sigmaA[i]^2 + sigmaL^2)
  G[i] = D[i] - E[i]
  G[i] = G[i]*(-0.5)*(1/(sigmaL^2 * sigmaA[i]^2))
  H[i] = exp(G[i])
  numerator[i] = H[i]/(sqrt(sigmaA[i]^2 + sigmaL^2))
  weightsQc[i] = (weights[i] * numerator[i])/(denominator)
}
muQc = c(0,0,0)
sigmaQc = muQc
for (i in 1:3){
  muQc[i] = (mu[i]*sigmaA[i]^2 + x*sigmaL^2)/(sigmaL^2 + sigmaA[i]^2)
  sigmaQc[i] = (sigmaL*sigmaA[i])/(sqrt(sigmaL^2 + sigmaA[i]^2))
}
plot(theta,densityValues(theta,muQc,sigmaQc,weightsQc), xlab = "RI", ylab = expression(pi(theta|x)),main = "Posterior distribution")
credibilityInterval = c(0,0)
for (i in 1:3){
  credibilityInterval[1] = credibilityInterval[1] + weightsQc[i]*(muQc[i] - 1.96*sqrt(sigmaQc[i]))
  credibilityInterval[2] = credibilityInterval[2] + weightsQc[i]*(muQc[i] + 1.96*sqrt(sigmaQc[i]))
}


# Q2 - d ------------------------------------------------------------------

#exponential term 
numeratorD = c(0,0,0)
denominatorD = c(0,0,0)
term = c(0,0,0)
weightsQd = weightsQc
muQd = c(0,0,0)
sigmaQd = muQd
for (i in 1:3){
  muQd[i] = (mu[i]*sigmaA[i]^2 + x*sigmaL^2)/(sigmaL^2 + sigmaA[i]^2)
  sigmaQd[i] = sqrt((2*sigmaL^2 * sigmaA[i]^2 + sigmaL^4)/(sigmaL^2 + sigmaA[i]^2))
}
plot(theta,densityValues(theta,muQd,sigmaQd,weightsQd), xlab = "RI", ylab = expression(pi(xnew|x)),main = "Posterior predictive density")


# Q2 - f ------------------------------------------------------------------

denominator = 0
x = theta
for (i in 1:3){
  A = (weights[i])/sqrt((sigmaA[i])^2 + sigmaL^2)
  B = exp((-0.5*((x-mu[i])^2))/((sigmaA[i])^2 + sigmaL^2))
  C = A*B
  denominator = denominator + C
}
{numerator = c(0,0,0)
  D = c(0,0,0)
  E = D
  G = E
  H = G
  weightsQf = numerator}
for (i in 1:3){
  D[i] = x^2 * sigmaL^2 + mu[i]^2*sigmaA[i]^2
  E[i] = (x * sigmaL^2 + mu[i]*sigmaA[i]^2)^2
  E[i] = E[i]/(sigmaA[i]^2 + sigmaL^2)
  G[i] = D[i] - E[i]
  G[i] = G[i]*(-0.5)*(1/(sigmaL^2 * sigmaA[i]^2))
  H[i] = exp(G[i])
  numerator[i] = H[i]/(sqrt(sigmaA[i]^2 + sigmaL^2))
  weightsQf[i] = (weights[i] * numerator[i])/(denominator)
}
muQc = c(0,0,0)
sigmaQc = muQc
for (i in 1:3){
  muQc[i] = (mu[i]*sigmaA[i]^2 + x*sigmaL^2)/(sigmaL^2 + sigmaA[i]^2)
  sigmaQc[i] = (sigmaL*sigmaA[i])/(sqrt(sigmaL^2 + sigmaA[i]^2))
}
plot(theta,densityValues(theta,muQc,sigmaQc,weightsQf), xlab = "RI", ylab = expression(pi(theta|x)),main = "Posterior distribution")
credibilityInterval = c(0,0)
for (i in 1:3){
  credibilityInterval[1] = credibilityInterval[1] + weightsQc[i]*(muQc[i] - 1.96*sqrt(sigmaQc[i]))
  credibilityInterval[2] = credibilityInterval[2] + weightsQc[i]*(muQc[i] + 1.96*sqrt(sigmaQc[i]))
}

theta2 <- seq(1.515,1.521,length = 100)
densityF <- rep(0,length(theta))
for (j in 1:length(theta)){
  denominator = 0
  x = theta2[j]
  for (i in 1:3){
    A = (weights[i])/sqrt((sigmaA[i])^2 + sigmaL^2)
    B = exp((-0.5*((x-mu[i])^2))/((sigmaA[i])^2 + sigmaL^2))
    C = A*B
    denominator = denominator + C
  }
  {numerator = c(0,0,0)
    D = c(0,0,0)
    E = D
    G = E
    H = G
    weightsQf = numerator}
  for (i in 1:3){
    D[i] = x^2 * sigmaL^2 + mu[i]^2*sigmaA[i]^2
    E[i] = (x * sigmaL^2 + mu[i]*sigmaA[i]^2)^2
    E[i] = E[i]/(sigmaA[i]^2 + sigmaL^2)
    G[i] = D[i] - E[i]
    G[i] = G[i]*(-0.5)*(1/(sigmaL^2 * sigmaA[i]^2))
    H[i] = exp(G[i])
    numerator[i] = H[i]/(sqrt(sigmaA[i]^2 + sigmaL^2))
    weightsQf[i] = (weights[i] * numerator[i])/(denominator)
  }
  muQf = c(0,0,0)
  sigmaQf = muQc
  for (i in 1:3){
    muQf[i] = (mu[i]*sigmaA[i]^2 + x*sigmaL^2)/(sigmaL^2 + sigmaA[i]^2)
    sigmaQf[i] = (sigmaL*sigmaA[i])/(sqrt(sigmaL^2 + sigmaA[i]^2))
  }
  densityF[j] <- densityValues(x,muQf,sigmaQf,weightsQf)
  densityF[j] <- 404.82/densityF[j]
}
plot(theta2,densityF,xlab = "Xc", ylab = "LR", main = "Likelihood ratio plot for Xc")
