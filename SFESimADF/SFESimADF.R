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

# set seed


# define sample size
iSample = 100
alpha   = c(0.9, 1)
beta    = c(-0.99, -0.90, 0, 0.90, 0.99)
pvec    = as.matrix(c(2,7))
x = rep(0, iSample)
set.seed(42)
epsilon = rnorm(iSample)
Simprocesses = function(alpha, beta){
  
  for(iCounter in 2:iSample){
    x[iCounter] = alpha * x[iCounter - 1] + beta * epsilon[iCounter - 1] + epsilon[iCounter]}
  return(x)}
params = cbind(sort(rep(alpha, length(beta))), rep(beta, length(alpha)))
processes       = mapply(Simprocesses, params[,1], params[,2])
# define function to simulate data

#


# apply adf test



adfforvaryingx = function(X){ADFtest(process, p)}
adfforvaryingp = function(p){ADFtest(process, p)}

Pvalues = matrix(NaN, nrow = length(pvec), ncol = ncol(processes))

ADFtest = function(x, p){
  pvalue = adf.test(x, alternative = c("stationary"), k = p)$p.value
  return(pvalue)}

adfforvaryingx = function(x){
  adfforvaryingp = function(p){
    ADFtest(x, p)}
  pvalues  = apply(pvec, 1, adfforvaryingp)
  return(pvalues)
}

colnames(processes) = c("StationaryProcess1", "StationaryProcess2", "StationaryProcess3", "StationaryProcess4", "StationaryProcess5",
                        "ExplosiveProcess1", "ExplosiveProcess2", "ExplosiveProcess3", "ExplosiveProcess4", "ExplosiveProcess5")
rownames(pvec) = c("2", "7")
pvaluesadf = apply(processes, 2, adfforvaryingx)
