# Solution_Q1.R
# This is the Solution file for Homework 6 (Question 1) of course:
# Bayesian Methods for Data Science (DATS 6450 - 11, Fall 2017)
# Data Science @ George Washington University
# Author: Yuxiao Huang

# Reference:
# Some of the code is from the book by Professor John K. Kruschke
# Please find the reference to and website of the book below:
# Kruschke, J. K. (2014). Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier
# https://sites.google.com/site/doingbayesiandataanalysis/


graphics.off()
rm(list=ls(all=TRUE))
source("DBDA2E-utilities.R")
require(rjags)
fileNameRoot="glucose" # for output filenames

#------------------------------------------------------------------------------- 
# Read the data 
data = read.csv("glucoseData.csv")

#------------------------------------------------------------------------------
# THE DATA.
yC = data$groupC
  
# Specify the data in a list, for later shipment to JAGS:
dataList = list(
  yC = yC,
  zC = sum(yC)
)
  
#-----------------------------------------------------------------------------
# THE MODEL.
modelString = "
model {
  zC ~ dbin( theta , length(yC) )

  theta ~ dbeta( omega[m]*(kappa[m]-2)+1 , (1-omega[m])*(kappa[m]-2)+1 ) 
  
  omega[1] = 0.5 # omega = (a - 1) / (a + b - 2)
  omegaTwo ~ dbeta( 50 , 50 )
  omega[2] = omegaTwo

  kappa[1] = 18 # kappa = a + b - 2
  kappaMinusTwo ~ dgamma( 0.01 , 0.01 )
  kappa[2] = kappaMinusTwo+ 2

  m ~ dcat( mPriorProb[] )
  mPriorProb[1] = .5
  mPriorProb[2] = .5
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

#-----------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Initial values of MCMC chains based on data:
# initsList = function() {
#   resampledZC = rbinom(1, size=length(yC) , prob=sum(yC)/length(yC) )
#   
#   thetaInit = resampledZC/length(yC)
#   
#   omegaTwoInit = thetaInit
#   
#   kappaMinusTwoInit = 100
#   
#   mInit = 2
#   
#   return( list( theta=thetaInit,
#                 omegaTwo = omegaTwoInit,
#                 kappaMinusTwo = kappaMinusTwoInit,
#                 m = mInit
#                 ))
# }

#-----------------------------------------------------------------------------
# RUN THE CHAINS
parameters = c("theta", "omega[2]", "kappa[2]", "m") 
adaptSteps = 1000             # Number of steps to "tune" the samplers.
burnInSteps = 1000           # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
numSavedSteps=50000          # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , # inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nPerChain , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcMat = as.matrix( codaSamples , chains=TRUE )
m = mcmcMat[,"m"]
theta = mcmcMat[,"theta"]

# Compute the proportion of m at each index value:
pM1 = sum( m == 1 ) / length( m )
pM2 = 1 - pM1

# Extract theta values for each model index:
thetaM1 = theta[ m == 1 ]
thetaM2 = theta[ m == 2 ]

# Plot histograms of sampled theta values for each model,
# with pM displayed.
openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,1,2,3),nrow=2,byrow=FALSE) , widths=c(1,2) )
plotPost( m , breaks=seq(0.9,2.1,0.2) , cenTend="mean" , xlab="m" , main="Model Index" )
plotPost( thetaM1 , 
          main=bquote( theta*" when m=1" * " ; p(m=1|D)" == .(signif(pM1,3)) ) , 
          cex.main=1.75 , xlab=bquote(theta) , xlim=c(0,1) )
plotPost( thetaM2 , 
          main=bquote( theta*" when m=2" * " ; p(m=2|D)" == .(signif(pM2,3)) ) , 
          cex.main=1.75 , xlab=bquote(theta) , xlim=c(0,1) )
saveGraph( file=paste0(fileNameRoot,"Post") , type="pdf" )