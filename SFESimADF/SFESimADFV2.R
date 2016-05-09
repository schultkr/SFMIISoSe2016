# remove variables
rm(list = ls())

# reset graphics
graphics.off()

# Install packages if not installed
libraries = c("tseries")
lapply(libraries, function(Samples) if (!(Samples %in% installed.packages())) {
  install.packages(Samples)
})

# Load packages
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# define number of simulations for ADF test
iSimulations = 1000
# define sample size
iSample      = 100
# define beta parmaters
beta = c(-0.99, -0.9, 0, 0.9, 0.99)
# define stationary models
statparams = cbind(rep(0.9, length(beta)), beta = c(-0.99, -0.9, 0, 0.9, 0.99))
# define explosive models
expparams  = cbind(rep(1, length(beta)), beta = c(-0.99, -0.9, 0, 0.9, 0.99))
#define lag vector
pvec         = as.matrix(3:11)
# define significance level
isiglevel    = 0.05
# define function to create simulated process
generateprocess = function(alpha, beta, iSample){
  epsilon = rnorm(iSample)
  x       = rep(0, iSample)
  for(iCounter in 2:iSample){
    x[iCounter] = alpha * x[iCounter - 1] + beta * epsilon[iCounter - 1] + epsilon[iCounter]}
  return(x)}

# adftest for varying p
adftestvaryp = function(pvec, x){
  adftestonep = function(p){return(adf.test(x, alternative = c("stationary"), k = p)$p.value)}
  pvalues     = apply(pvec, 1, adftestonep)
  rejection   = as.numeric(pvalues < isiglevel)
  return(rejection)
}

# define function to simulate rejection probabilities
ADFSimtest = function(alpha, beta){
  rejection = matrix(rep(NaN, iSimulations, length(pvec)), nrow = iSimulations, ncol = length(pvec))
  for(iCounterSim in 1:iSimulations){
    x                        = generateprocess(alpha, beta, iSample)
    rejection[iCounterSim, ] = adftestvaryp(pvec, x)}
  rejectionprob = colMeans(rejection)
  return(rejectionprob)
}

# simulate tests
teststat = mapply(ADFSimtest, statparams[, 1], statparams[, 2])
testexp  = mapply(ADFSimtest, expparams[, 1], expparams[, 2])

# plot results
matplot(teststat, type = "l", lwd = 3, ylab = "Rejection Probability", xlab = "Lags", main = "Level of ADF")
matplot(testexp, type = "l", lwd = 3, ylab = "Rejection Probability", xlab = "Lags", main = "Power of ADF")

# create latex table
tablestationary = round(teststat[c(1, 9), ], 2)
tableexplosive  = round(testexp[c(1, 9), ], 2)
colnames(tablestationary) = sapply(beta, as.character)
rownames(tablestationary) = sapply(c(3, 11), as.character)
