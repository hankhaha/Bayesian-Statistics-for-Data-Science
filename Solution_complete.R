# Solution.R
# This is the Solution file for Exercise 4 of course:
# Bayesian Methods for Data Science (DATS 6450 - 11, Spring 2018)
# Data Science @ George Washington University
# Author: Yuxiao Huang

# Reference:
# Some of the code is from the book by Professor John K. Kruschke
# Please find the reference to and website of the book below:
# Kruschke, J. K. (2014). Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier
# https://sites.google.com/site/doingbayesiandataanalysis/

# Students should implement the following function in Solution.R:
# genMCMC: generates the MCMC chain

source("DBDA2E-utilities.R")

#===============================================================================
genMCMC = function( data, numSavedSteps=50000 , saveName=NULL , thinSteps=1 ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  require(rjags)
  require(runjags)
  
  #-----------------------------------------------------------------------------
  # THE DATA.
  yA = data$groupA
  yB = data$groupB
  yC = data$groupC
  
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    yA = yA,
    yB = yB,
    yC = yC,
    zA = sum(yA),
    zB = sum(yB),
    zC = sum(yC)
  )
  
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    zA ~ dbin( thetaA , length(yA) )
    zB ~ dbin( thetaB , length(yB) )
    zC ~ dbin( thetaC , length(yC) )

    thetaA ~ dbeta( omega*(kappa-2)+1 , (1-omega)*(kappa-2)+1 )
    thetaB ~ dbeta( omega*(kappa-2)+1 , (1-omega)*(kappa-2)+1 )
    thetaC ~ dbeta( 10 , 10 )

    omega ~ dbeta( 10 , 10 )

    # kappaMinusTwo ~ dgamma( 0.01 , 0.01 )
    kappaMinusTwo ~ dgamma( 1000000 , 1000 )

    kappa = kappaMinusTwo + 2
  }
  " # close quote for modelString
  writeLines( modelString , con="TEMPmodel.txt" )
  
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  initsList = function() {
    resampledZA = rbinom(1, size=length(yA) , prob=sum(yA)/length(yA) )
    resampledZB = rbinom(1, size=length(yB) , prob=sum(yB)/length(yB) )
    resampledZC = rbinom(1, size=length(yC) , prob=sum(yC)/length(yC) )
    
    thetaInitA = resampledZA/length(yA)
    thetaInitB = resampledZB/length(yB)
    thetaInitC = resampledZC/length(yC)
    
    omegaInit = mean(c(thetaInitA, thetaInitB))
    
    kappaMinusTwoInit = 100
    
    return( list( thetaA=thetaInitA,
                  thetaB=thetaInitB,
                  thetaC=thetaInitC,
                  
                  omega = omegaInit,
                
                  kappaMinusTwo = kappaMinusTwoInit
                  ))
  }
  
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "thetaA", "thetaB", "thetaC", "omega", "kappa") 
  adaptSteps = 500             # Number of steps to adapt the samplers
  burnInSteps = 500            # Number of steps to burn-in the chains
  
  useRunjags = TRUE
  if ( useRunjags ) {
    runJagsOut <- run.jags( method=runjagsMethod ,
                            model="TEMPmodel.txt" , 
                            monitor=parameters , 
                            data=dataList ,  
                            inits=initsList , 
                            n.chains=nChains ,
                            adapt=adaptSteps ,
                            burnin=burnInSteps , 
                            sample=ceiling(numSavedSteps/nChains) ,
                            thin=thinSteps ,
                            summarise=FALSE ,
                            plots=FALSE )
    codaSamples = as.mcmc.list( runJagsOut )
  } else {
    # Create, initialize, and adapt the model:
    jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                            n.chains=nChains , n.adapt=adaptSteps )
    # Burn-in:
    cat( "Burning in the MCMC chain...\n" )
    update( jagsModel , n.iter=burnInSteps )
    # The saved MCMC chain:
    cat( "Sampling final MCMC chain...\n" )
    codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                                n.iter=ceiling(numSavedSteps*thinSteps/nChains), 
                                thin=thinSteps )
  }  
  
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
}

#===============================================================================
plotMCMC = function( codaSamples , data ,
                     compVal=0.5 , rope=NULL , 
                     diffSList=NULL , diffCList=NULL , 
                     compValDiff=0.0 , ropeDiff=NULL , 
                     saveName=NULL , saveType="jpg" ) {
  # Now plot the posterior:
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  
  # kappa:
  parNames = sort(grep("kappa",colnames(mcmcMat),value=TRUE))
  nPanels = length(parNames)
  nCols = 1
  nRows = ceiling(nPanels/nCols)
  openGraph(width=2.5*nCols,height=2.0*nRows)
  par( mfcol=c(nRows,nCols) )
  par( mar=c(3.5,1,3.5,1) , mgp=c(2.0,0.7,0) )
  #xLim = range( mcmcMat[,parNames] )
  xLim=quantile(mcmcMat[,parNames],probs=c(0.000,0.995))
  #mainLab = c(levels(myData[[cName]]),"Overall")
  mainIdx = 0
  for ( parName in parNames ) {
    mainIdx = mainIdx+1
    postInfo = plotPost( mcmcMat[,parName] , compVal=compVal , ROPE=rope ,
                         xlab=bquote(.(parName)) , cex.lab=1.25 , 
                         main=NULL , cex.main=1.5 ,
                         xlim=xLim , border="skyblue" )
  }  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Kappa",sep=""), type=saveType)
  }
  
  # omega:
  parNames = sort(grep("omega",colnames(mcmcMat),value=TRUE))
  nPanels = length(parNames)
  nCols = 4
  nRows = ceiling(nPanels/nCols)
  openGraph(width=2.5*nCols,height=2.0*nRows)
  par( mfcol=c(nRows,nCols) )
  par( mar=c(3.5,1,3.5,1) , mgp=c(2.0,0.7,0) )
  #xLim = range( mcmcMat[,parNames] )
  xLim=quantile(mcmcMat[,parNames],probs=c(0.001,0.999))
  #mainLab = c(levels(myData[[cName]]),"Overall")
  mainIdx = 0
  for ( parName in parNames ) {
    mainIdx = mainIdx+1
    postInfo = plotPost( mcmcMat[,parName] , compVal=compVal , ROPE=rope ,
                         xlab=bquote(.(parName)) , cex.lab=1.25 , 
                         main=NULL , cex.main=1.5 ,
                         xlim=xLim , border="skyblue" )
  }  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Omega",sep=""), type=saveType)
  }
  
  # theta:
  parNames = sort(grep("theta",colnames(mcmcMat),value=TRUE))
  nPanels = length(parNames)
  nCols = 4
  nRows = ceiling(nPanels/nCols)
  openGraph(width=2.5*nCols,height=2.0*nRows)
  par( mfcol=c(nRows,nCols) )
  par( mar=c(3.5,1,3.5,1) , mgp=c(2.0,0.7,0) )
  #xLim = range( mcmcMat[,parNames] )
  xLim=quantile(mcmcMat[,parNames],probs=c(0.001,0.999))
  #mainLab = c(levels(myData[[cName]]),"Overall")
  mainIdx = 0
  for ( parName in parNames ) {
    mainIdx = mainIdx+1
    postInfo = plotPost( mcmcMat[,parName] , compVal=compVal , ROPE=rope ,
                         xlab=bquote(.(parName)) , cex.lab=1.25 , 
                         main=NULL , cex.main=1.5 ,
                         xlim=xLim , border="skyblue" )
  }  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Omega",sep=""), type=saveType)
  }
  
  # Plot individual theta's and differences:
  if ( !is.null(diffSList) ) {
    for ( compIdx in 1:length(diffSList) ) {
      diffSVec = diffSList[[compIdx]]
      Nidx = length(diffSVec)
      openGraph(width=2.5*Nidx,height=2.0*Nidx)
      par( mfrow=c(Nidx,Nidx) )
      xLim = range(c( compVal, rope,
                      mcmcMat[,diffSVec] ))
      for ( t1Idx in 1:Nidx ) {
        for ( t2Idx in 1:Nidx ) {
          parName1 = diffSVec[t1Idx]
          parName2 = diffSVec[t2Idx]
          if ( t1Idx > t2Idx) {  
            # plot.new() # empty plot, advance to next
            par( mar=c(3,3,3,1) , mgp=c(2.0,0.7,0) , pty="s" )
            nToPlot = 700
            ptIdx = round(seq(1,chainLength,length=nToPlot))
            plot ( mcmcMat[ptIdx,parName2] , mcmcMat[ptIdx,parName1] , 
                   cex.main=1.25 , cex.lab=1.25 , 
                   xlab=diffSVec[t2Idx] , 
                   ylab=diffSVec[t1Idx] , 
                   col="skyblue" )
            abline(0,1,lty="dotted")
          } else if ( t1Idx == t2Idx ) {
            par( mar=c(3,1.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1] , 
                                 compVal=compVal , ROPE=rope , 
                                 cex.main=1.25 , cex.lab=1.25 , 
                                 xlab=bquote(.(parName1)) ,
                                 main=diffSVec[t1Idx] ,  
                                 xlim=xLim )
          } else if ( t1Idx < t2Idx ) {
            par( mar=c(3,1.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                                 compVal=compValDiff , ROPE=ropeDiff , 
                                 cex.main=1.25 , cex.lab=1.25 , 
                                 xlab=bquote("Difference of "*omega*"'s") , 
                                 main=paste( diffSVec[t1Idx] ,
                                   "-",diffSVec[t2Idx] ) )
          }
        }
      }
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,"OmegaDiff",compIdx), type=saveType)
      }
    }
  }
}