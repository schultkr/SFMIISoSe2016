# === preprocessing ===
# remove variables
rm(list = ls())

# reset graphics
graphics.off()
set.seed(48)
# Install packages if not installed
libraries = c("tseries")
lapply(libraries, function(Samples) if (!(Samples %in% installed.packages())) {
    install.packages(Samples)
})

# Load packages
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# === Input Parameters ===
# define number of simulations for ADF test
iSimulations = 1000

# define sample size
iSample = 100

# define alpha parameters
alpha = c(0.9, 1)

# define beta parmaters
beta = c(-0.99, -0.9, 0, 0.9, 0.99)

# define stationary models
statparams = cbind(rep(alpha[1], length(beta)), beta = c(-0.99, -0.9, 0, 0.9, 0.99))

# define explosive models
expparams = cbind(rep(alpha[2], length(beta)), beta = c(-0.99, -0.9, 0, 0.9, 0.99))

# define lag vector
pvec = as.matrix(3:11)

# define significance level
isiglevel = 0.05

# === Define Functions ===
# define function to create simulated process
generateprocess = function(alpha, beta, iSample) {
    epsilon = rnorm(iSample)
    x       = rep(0, iSample)
    for (iCounter in 2:iSample) {
        x[iCounter] = alpha * x[iCounter - 1] + beta * epsilon[iCounter - 1] + epsilon[iCounter]
    }
    return(x)
}

# adftest for varying p
adftestvaryp = function(pvec, x) {
    adftestonep = function(p) {
        return(adf.test(x, alternative = c("stationary"), k = p)$p.value)
    }
    pvalues     = apply(pvec, 1, adftestonep)
    rejection   = as.numeric(pvalues < isiglevel)
    return(rejection)
}

# define function to simulate rejection probabilities
ADFSimtest = function(alpha, beta) {
    res = matrix(rep(NaN, iSimulations, length(pvec)), nrow = iSimulations, ncol = length(pvec))
    for (iSim in 1:iSimulations) {
        x           = generateprocess(alpha, beta, iSample)
        res[iSim, ] = adftestvaryp(pvec, x)
    }
    rejectionprob = colMeans(res)
    return(rejectionprob)
}

# === Main Computation ===
# simulate tests
teststat = mapply(ADFSimtest, statparams[, 1], statparams[, 2])
testexp  = mapply(ADFSimtest, expparams[, 1], expparams[, 2])


# plot power of test
par(mfrow = c(1, 2), cex.lab = 1.1)
matplot(pvec, teststat * 100, type = "l", lwd = 3, ylab = "Rejection Probability", xlab = "Lags", 
        main = "Power of ADF Test", col = c("black", "red3", "blue3", "green3", "magenta3"), 
        xlim = c(min(pvec), max(pvec)), ylim = c(0, 100))

# plot level of test
matplot(pvec, testexp * 100, type = "l", lwd = 3, ylab = "Rejection Probability", xlab = "Lags", 
        main = "Level of ADF Test", col = c("black", "red3", "blue3", "green3", "magenta3"), 
        xlim = c(min(pvec), max(pvec)), ylim = c(0, 100))

# round results and use only p = 3 and p = 11
tablestationary = round(teststat[c(1, length(pvec)), ], digits = 3)
tableexplosive  = round(testexp[c(1, length(pvec)), ], digits = 3)

tablehelp       = rbind(c(" ", " ", "beta", " ", " "), sapply(beta, as.character), tablestationary)
tablestatprint  = cbind(c(" ", "alpha", alpha[1], " "), c(" ", "p", pvec[1], pvec[length(pvec)]), tablehelp)

tablehelp       = rbind(c(" ", " ", "beta", " ", " "), sapply(beta, as.character), tableexplosive)
tableexpprint   = cbind(c(" ", "alpha", alpha[2], " "), c(" ", "p", pvec[1], pvec[length(pvec)]), tablehelp)

# print tables
options(digits = 3)
cat("\014")
# table for power of test
tablestatprint
# table for level of test
tableexpprint

