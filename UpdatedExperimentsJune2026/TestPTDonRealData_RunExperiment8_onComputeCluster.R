#setwd("~/")

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
SoftwareImplementationPath <- paste0(getwd(),"/PTD_Boot_Implementation/")
source(paste0(SoftwareImplementationPath,"PTDBootModularized.R"))
source(paste0(SoftwareImplementationPath,"BootstrapAndCalcEsts.R"))
source(paste0(SoftwareImplementationPath,"EstBetaGammaCalib_andVCOV.glm.R"))
source(paste0(SoftwareImplementationPath,"FormatAndFitReg.R"))


library(doParallel)
runSimsRealData <- function(FormattedDatInp2,N_sampInp,nLabelTargInp,nSims=100,regTypeInp="linear",outcomeVarName="Y",BootstrapSchemesInp=c("FullBootPercentile","QuickConvolution","CLTBased"),tuningSchemesInp=c("Optimal","Diagonal","None"),alphaInp=0.1,clusteredLabelling=F,useClusterIDs=F,tauForQuantileReg=NULL,nBootInferenceInp2 = 2000,ncores=35){
  
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
    
    suppressWarnings(PTDBootModularized(ProxyDat = FormattedDatSubsamp$ProxyDat,GoodDat = CurrLabelling$GoodDatWMissing,PiLabel = CurrLabelling$piLabelRescaled,clusterID=ClustersUse ,BootstrapScheme=BootstrapSchemesInp,nBootInference = nBootInferenceInp2,nBootTune = 0,TuningScheme = tuningSchemesInp,RegType = regTypeInp,OutcomeVarName = outcomeVarName,alpha = alphaInp,tauQuantReg = tauForQuantileReg)) 
    
  }
  totalTime <- (proc.time()[3]-startTime)/60
  file.remove("progress.log")
  stopImplicitCluster()
  
  return(list(SimResultList=SimResultList,N_sampInp=N_sampInp,nLabelTargInp=nLabelTargInp,RegTypeUsed=regTypeInp,alphaUsed=alphaInp,tauForQuantileRegUsed=tauForQuantileReg,totalTime=totalTime)) 
}



SimsResultPath <- paste0(getwd(),"/SimulationResults/")






#1500 simulations were run already. Running remaining 500 on computing cluster
load(paste0(getwd(),"/Datasets/CensusStratEducationFormatted.RData"))

CensusOrdinalRegSimsB2 <- runSimsRealData(FormattedDatInp = EducationDatFormatted,nSims = 500,N_sampInp = 6000,nLabelTargInp = 1000,regTypeInp = "ordinal",outcomeVarName = "Education",tuningSchemesInp = c("Optimal","Diagonal","None"),BootstrapSchemesInp = c("FullBootPercentile","QuickConvolution"))

save(CensusOrdinalRegSimsB2,file = paste0(SimsResultPath,"CensusOrdinalRegSimsBatch2.RData"))

