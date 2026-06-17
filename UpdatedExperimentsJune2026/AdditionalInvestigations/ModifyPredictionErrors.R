#setwd("~/AdditionalInvestigations")





RandomSubsetFixedSize <- function(FormattedDatInp,N_samp,clusterSample=F){
  
  FormattedDatSubset <- list()
  
  if(clusterSample){
    uniqueClusters <- unique(FormattedDatInp$ClusterIDs)
    NClusters <- length(uniqueClusters)
    ExpectedNClustSubsamp <- NClusters* N_samp/nrow(FormattedDatInp$GoodDat)
    rem <- ExpectedNClustSubsamp-floor(ExpectedNClustSubsamp)
    NClustSubsamp <- floor(ExpectedNClustSubsamp)+ sample(c(0,1),size = 1,prob = c(1-rem,rem))
    ClusterIDsSubsamp <- sample(uniqueClusters,size = NClustSubsamp,replace = F)
    idxSubsample <- which(FormattedDatInp$ClusterIDs %in% ClusterIDsSubsamp)
  } else{
    MSuper <- nrow(FormattedDatInp$GoodDat)
    idxSubsample <- sample(1:MSuper,size = N_samp,replace = F) #If you sample with replacement simulations should have more precise coverage
  }
  
  FormattedDatSubset$GoodDat <- FormattedDatInp$GoodDat[idxSubsample,]
  FormattedDatSubset$ProxyDat <- FormattedDatInp$ProxyDat[idxSubsample,]
  
  
  
  if(!is.null(FormattedDatInp$PiLabelUnsc)){FormattedDatSubset$PiLabelUnsc <- FormattedDatInp$PiLabelUnsc[idxSubsample]} 
  else{ FormattedDatSubset$PiLabelUnsc <- NULL}
  if(!is.null(FormattedDatInp$ClusterIDs)){FormattedDatSubset$ClusterIDs <- FormattedDatInp$ClusterIDs[idxSubsample]} 
  else{ FormattedDatSubset$ClusterIDs <- NULL}
  
  return(FormattedDatSubset)
}



LabelRandomSubset <- function(GoodDatInp,ProxyDatInp,piLabelInpUsnc=NULL,nLabelTarget){
  if(is.null(piLabelInpUsnc)){
    piLabelRescaled <- rep(nLabelTarget/nrow(GoodDatInp),nrow(GoodDatInp))
  } else {
    piLabelRescaled <- nLabelTarget*piLabelInpUsnc/sum(piLabelInpUsnc)
  }
  Idx2Label <- I(runif(n = nrow(GoodDatInp)) < piLabelRescaled)
  GoodDatWMissing <- GoodDatInp
  GoodDatWMissing[!Idx2Label,] <- NA
  return(list(GoodDatWMissing=GoodDatWMissing,piLabelRescaled=piLabelRescaled))
}


LabelRandomSubsetofClusters <- function(GoodDatInp,ProxyDatInp,piLabelInpUsnc=NULL,clustersInp,nLabelTarget){
  
  if(is.null(piLabelInpUsnc)){
    uniqueClusterIDs <- unique(clustersInp)
    piClusterRescaled <- rep(nLabelTarget/nrow(GoodDatInp),length(uniqueClusterIDs))
    piLabelRescaled <- rep(nLabelTarget/nrow(GoodDatInp),nrow(GoodDatInp))
  } else {
    clusterInp2 <- clustersInp
    piClusterDf <- data.frame(piLabelInpUsnc,clusterInp2) %>% group_by(clusterInp2) %>% summarise(piCluster=mean(piLabelInpUsnc),ClusterSize=n())
    piClusterRescaled <- nLabelTarget*piClusterDf$piCluster/sum(piClusterDf$piCluster*piClusterDf$ClusterSize)
    piLabelRescaled <- nLabelTarget*piLabelInpUsnc/sum(piLabelInpUsnc)
    uniqueClusterIDs <- piClusterDf[['clusterInp2']]
  }
  
  
  Clusters2Label <- I(runif(n = length(uniqueClusterIDs)) < piClusterRescaled)
  Idx2Label <- rep(NA,nrow(GoodDatInp))
  for(i in 1:length(uniqueClusterIDs)){
    Idx2Label[clustersInp==uniqueClusterIDs[i]] <- Clusters2Label[i]
  }
  GoodDatWMissing <- GoodDatInp
  GoodDatWMissing[!Idx2Label,] <- NA
  return(list(GoodDatWMissing=GoodDatWMissing,piLabelRescaled=piLabelRescaled))
}




library(dplyr)
SoftwareImplementationPath <- paste0(dirname(getwd()),"/PTD_Boot_Implementation/")
source(paste0(SoftwareImplementationPath,"PTDBootModularized.R"))
source(paste0(SoftwareImplementationPath,"BootstrapAndCalcEsts.R"))
source(paste0(SoftwareImplementationPath,"EstBetaGammaCalib_andVCOV.glm.R"))
source(paste0(SoftwareImplementationPath,"FormatAndFitReg.R"))

library(doParallel)
runSimsRealData <- function(FormattedDatInp2,N_sampInp,nLabelTargInp,nSims=100,regTypeInp="linear",outcomeVarName="Y",BootstrapSchemesInp=c("FullBootPercentile","QuickConvolution","CLTBased"),tuningSchemesInp=c("Optimal","Diagonal","None"),alphaInp=0.1,clusteredLabelling=F,useClusterIDs=F,tauForQuantileReg=NULL,sigmaNightlightNoise=0,ncores=78){
  
  
  
  registerDoParallel(ncores)
  
  
  writeLines("", "progress.log")
  list()
  currTime <- proc.time()[3]
  startTime <- currTime
  SimResultList <- foreach(i = 1:nSims) %dopar% {
    
    cat(paste0("Starting Simulation ", i, " out of ", nSims,
               " | ", format(Sys.time(), "%H:%M:%S"), "\n"),
        file = "progress.log", append = TRUE)
  

    FormattedDatSubsamp <- RandomSubsetFixedSize(FormattedDatInp = FormattedDatInp2,N_samp = N_sampInp,clusterSample = clusteredLabelling)
    
    if(clusteredLabelling){
      CurrLabelling <- LabelRandomSubsetofClusters(GoodDatInp = FormattedDatSubsamp$GoodDat,
                                                   ProxyDatInp = FormattedDatSubsamp$ProxyDat,
                                                   piLabelInpUsnc = FormattedDatSubsamp$PiLabelUnsc,
                                                   clustersInp = FormattedDatSubsamp$ClusterIDs,
                                                   nLabelTarget = nLabelTargInp)
    }else{
      CurrLabelling <- LabelRandomSubset(GoodDatInp = FormattedDatSubsamp$GoodDat,
                                         ProxyDatInp = FormattedDatSubsamp$ProxyDat,
                                         piLabelInpUsnc = FormattedDatSubsamp$PiLabelUnsc,
                                         nLabelTarget = nLabelTargInp)
    }
    
    if(useClusterIDs){
      ClustersUse <- FormattedDatSubsamp$ClusterIDs
    } else{
      ClustersUse <- NULL
    }
    
    if(sigmaNightlightNoise==0){
      proxyDatCurrSim <- FormattedDatSubsamp$ProxyDat
    } else { #Adding noise to nightlights
      proxyDatCurrSim <- FormattedDatSubsamp$ProxyDat
      NightlightsPredWithNoise <- FormattedDatSubsamp$ProxyDat$Nightlights+sigmaNightlightNoise*rnorm(n = nrow(proxyDatCurrSim))
      proxyDatCurrSim$Nightlights <- NightlightsPredWithNoise
    }
    
    suppressWarnings(PTDBootModularized(ProxyDat = proxyDatCurrSim,GoodDat = CurrLabelling$GoodDatWMissing,PiLabel = CurrLabelling$piLabelRescaled,clusterID=ClustersUse ,BootstrapScheme=BootstrapSchemesInp,nBootInference = 2000,nBootTune = 0,TuningScheme = tuningSchemesInp,RegType = regTypeInp,OutcomeVarName = outcomeVarName,alpha = alphaInp,tauQuantReg = tauForQuantileReg)) 
    
  }
  totalTime <- (proc.time()[3]-startTime)/60
  file.remove("progress.log")
  stopImplicitCluster()
  
  return(list(SimResultList=SimResultList,N_sampInp=N_sampInp,nLabelTargInp=nLabelTargInp,RegTypeUsed=regTypeInp,alphaUsed=alphaInp,tauForQuantileRegUsed=tauForQuantileReg,totalTime=totalTime)) 
}




ExtractSummaryTableSingleExp <- function(SimResultsInp,CoefNamesInp,ErrorParam,MethodNameMat,TrueBeta,TrueGamma){
  
  MethodNameVec <-  MethodNameMat[,1] 
  MethodNameFormatted <- MethodNameMat[,2] 
  df1TuneMethodList <- list()
  for(methId in 1:length(MethodNameVec)){
    nCoeff <- nrow(CoefNamesInp)
    CIWidthsTabCurr <- matrix(NA,ncol = length(SimResultsInp),nrow = nCoeff)
    CoversMat <- matrix(NA,ncol = length(SimResultsInp),nrow = nCoeff)
    
    for(j in 1:length(SimResultsInp)){
      CIsCurrSim <- SimResultsInp[[j]]$CIList[[MethodNameVec[methId]]]
      if((sum(rownames(CIsCurrSim) != CoefNamesInp[,1])==0) | (sum(rownames(CIsCurrSim)  != names(TrueBeta))==0) ){
        CIWidthsTabCurr[,j] <- CIsCurrSim[,'CI_ub']-CIsCurrSim[,'CI_lb']
        CoversMat[,j] <- I((TrueBeta<CIsCurrSim[,'CI_ub']) & (TrueBeta > CIsCurrSim[,'CI_lb']))
      } else{print("Warning:Coefficient name misalignment")}
    }
    
    CIWidths_mean <- rowMeans(CIWidthsTabCurr)
    CIWidths_sd <- rep(NA,nCoeff)
    for(i in 1:nCoeff){CIWidths_sd[i] <- sd(CIWidthsTabCurr[i,])}
    EmpCoverage <- rowMeans(CoversMat)
    
    dfOutCurr <- data.frame(rep(ErrorParam,nCoeff),CoefNamesInp[,2],rep(MethodNameFormatted[methId],nCoeff),
                            CIWidths_mean,CIWidths_sd,EmpCoverage,TrueBeta,TrueGamma)
    names(dfOutCurr) <- c("ErrorParameter","Coefficient","Method","CIWidths_mean","CIWidths_sd","Coverage","TrueBeta","TrueGamma")
    df1TuneMethodList[[methId]] <- dfOutCurr
  }
  dfOut <- df1TuneMethodList[[1]]
  for(methId in 2:length(MethodNameVec)){dfOut <- rbind(dfOut,df1TuneMethodList[[methId]])}
  return(dfOut)
}


load(paste0(dirname(getwd()),"/Datasets/","HousingPriceFormatted.RData"))

#Extract info used for summary statistic table
TrueBetaHousingPriceLM <- lm(formula = HousingPrice ~ .,data = HousingPriceFormatted$GoodDat)$coefficients
TrueGammaHousingPriceLM <- lm(formula = HousingPrice ~ .,data = HousingPriceFormatted$ProxyDat)$coefficients
CoefNameMatHousingPrice <- cbind(c("(Intercept)","Income","Nightlights","RoadLength"),
                                 c("Intercept","Income","Nightlights","Road Length"))


#Two types of modification
# 1) \tilde{Y}'=Y+t*(\tilde{Y}-Y) #Prediction error in same direction but larger

RescaleNightlightError <- function(dfInp,tInp){
  NightlightsGood <- dfInp$GoodDat$Nightlights
  NightlightsProxyOrig <- dfInp$ProxyDat$Nightlights
  dfOut <- dfInp
  dfOut$ProxyDat$Nightlights <- NightlightsGood + tInp*(NightlightsProxyOrig-NightlightsGood)
  return(dfOut)
}

# 2) \tilde{Y}'=Y+t*(\tilde{Y}-Y) #This is implemented by varying sigmaNightlightNoise parameter in the runSimsRealData function



#Storing names of outputs considered
MethodNameMatConvBootBased <- cbind(c(paste0("QuickConvolution_",c("Optimal","Diagonal","None")),"classical"),
                                    c("PTD (Optimal Tuning)","PTD (Optimal Diagonal Tuning)","PTD (No Tuning)","Classical"))

SimsResultPath <- paste0(getwd(),"/Results/IncreaseError/")
nSimsPerExp <- 2000

#Running Extrapolated (or interpolated error) experiments

tVec <- c(0,0.25,0.5,0.75,seq(1,6,by = 0.5))

ModifyNightlightProxyErrorInterpAndExtrap <- list()

for(t in 1:length(tVec)){
  print(paste0("t=",tVec[t]))
  FormattedDatModif <- RescaleNightlightError(dfInp = HousingPriceFormatted,tInp = tVec[t])
  HousingPriceSimsCurr <- runSimsRealData(FormattedDatInp = FormattedDatModif,nSims = nSimsPerExp,N_sampInp = 5000,nLabelTargInp = 500,regTypeInp = "linear",outcomeVarName = "HousingPrice",BootstrapSchemesInp = c("QuickConvolution"),sigmaNightlightNoise = 0)
  
  
  ModifyNightlightProxyErrorInterpAndExtrap[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = HousingPriceSimsCurr$SimResultList,CoefNamesInp = CoefNameMatHousingPrice,ErrorParam  = tVec[t],TrueBeta = TrueBetaHousingPriceLM ,TrueGamma = TrueGammaHousingPriceLM,MethodNameMat = MethodNameMatConvBootBased)
}


save(tVec,nSimsPerExp,ModifyNightlightProxyErrorInterpAndExtrap, file = paste0(SimsResultPath,"/ExrapolateInterpolateErrorQuickBoot.RData"))

## Running adding noise experiments

addNightlightNoiseSims <- list()


sigmaVec <- seq(0,5,by = 0.5)

for(sigIdx in 1:length(sigmaVec)){
  print(paste0("sigma= ",sigmaVec[sigIdx]))
  HousingPriceSimsCurr <- runSimsRealData(FormattedDatInp = HousingPriceFormatted,nSims = nSimsPerExp,N_sampInp = 5000,nLabelTargInp = 500,regTypeInp = "linear",outcomeVarName = "HousingPrice",BootstrapSchemesInp = c("QuickConvolution"),sigmaNightlightNoise = sigmaVec[sigIdx],freqToc = 50)
  
  addNightlightNoiseSims[[sigIdx]] <- ExtractSummaryTableSingleExp(SimResultsInp = HousingPriceSimsCurr$SimResultList,CoefNamesInp = CoefNameMatHousingPrice,ErrorParam  = sigmaVec[sigIdx],TrueBeta = TrueBetaHousingPriceLM ,TrueGamma = TrueGammaHousingPriceLM,MethodNameMat = MethodNameMatConvBootBased)
}

save(sigmaVec,nSimsPerExp,addNightlightNoiseSims, file = paste0(SimsResultPath,"/IncreaseNoiseSimsQuickBoot.RData"))
