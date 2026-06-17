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
runSimsRealData <- function(FormattedDatInp2,N_sampInp,nLabelTargInp,nSims=100,regTypeInp="linear",outcomeVarName="Y",BootstrapSchemesInp=c("FullBootPercentile","QuickConvolution","CLTBased"),tuningSchemesInp=c("Optimal","Diagonal","None"),alphaInp=0.1,clusteredLabelling=F,useClusterIDs=F,tauForQuantileReg=NULL,nBootInferenceInp2 = 2000,ncores=70){
  
  registerDoParallel(ncores)
  
  
  writeLines("", "progress.log")
  list()
  currTime <- proc.time()[3]
  startTime <- currTime
  SimResultList <- foreach(i = 1:nSims) %dopar% {
    
    cat(paste0("Starting Simulation ", i, " out of ", nSims, "sim with n=", nLabelTargInp,
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


library(dplyr)
ExtractSummaryTableSingleExp <- function(SimResultsInp,CoefNamesInp,nInp,MethodNameMat,TrueBeta){
  
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
    
    CIWidths_median <- rep(NA,nCoeff)
    CIWidths_sd <- rep(NA,nCoeff)
    for(i in 1:nCoeff){
      CIWidths_median[i] <- median(CIWidthsTabCurr[i,])
      CIWidths_sd[i] <- sd(CIWidthsTabCurr[i,])
    }
    EmpCoverage <- rowMeans(CoversMat)
    
    
    dfOutCurr <- data.frame(rep(nInp,nCoeff),CoefNamesInp[,2],rep(MethodNameFormatted[methId],nCoeff),
                            CIWidths_mean,CIWidths_median,CIWidths_sd,EmpCoverage,TrueBeta)
    names(dfOutCurr) <- c("n","Coefficient","Method","CIWidths_mean","CIWidths_median","CIWidths_sd","Coverage","TrueBeta")
    df1TuneMethodList[[methId]] <- dfOutCurr
  }
  dfOut <- df1TuneMethodList[[1]]
  for(methId in 2:length(MethodNameVec)){dfOut <- rbind(dfOut,df1TuneMethodList[[methId]])}
  return(dfOut)
}


#Load datasets
load(paste0(dirname(getwd()),"/Datasets/","HousingPriceFormatted.RData"))
load(paste0(dirname(getwd()),"/Datasets/","TreecoverBinarizedFormatted.RData"))
load(paste0(dirname(getwd()),"/Datasets/","TreecoverFormatted.RData"))
load(paste0(dirname(getwd()),"/Datasets/","AlphaFoldFormatted.RData"))



MethodNameMatPTD <- cbind(c("FullBootPercentile_Optimal","FullBootPercentile_Diagonal","FullBootPercentile_None","QuickConvolution_Optimal","QuickConvolution_Diagonal","QuickConvolution_None","classical"),c("PTD (Optimal Tuning)","PTD (Optimal Diagonal Tuning)","PTD (No Tuning)","PTD Speedup (Optimal Tuning)","PTD Speedup (Optimal Diagonal Tuning)","PTD Speedup (No Tuning)","Classical"))


CoefNameMatAlphaFold <- cbind(c("(Intercept)","ubiquitinated","acetylated","`ubiq:acet_interaction`"),
         c("Intercept","Ubiquitinated","Acetylated","Interaction"))

CoefNameMatTreecover <- cbind(c("(Intercept)","Elevation","Population"),
                              c("Intercept","Elevation","Population"))


CoefNameMatHousingPrice <- cbind(c("(Intercept)","Income","Nightlights","RoadLength"),
                                 c("Intercept","Income","Nightlights","RoadLength"))




#Calculating True Regression Coefficients

TrueBetaAlphaFoldLogistic <- glm(formula = IDR ~ .,data = AlphaFoldFormatted$GoodDat, family =  binomial(link = "logit"))$coefficients

TrueBetaHousingPriceLM <- lm(formula = HousingPrice ~ .,data = HousingPriceFormatted$GoodDat)$coefficients

library(quantreg)
TrueBetaHousingPriceQR <-  rq(formula = HousingPrice ~ .,data = HousingPriceFormatted$GoodDat,tau = 0.5)$coefficients


TrueBetaTreecoverLM <- lm(formula = Treecover ~ .,data = TreecoverFormatted$GoodDat)$coefficients

TrueBetaForestedLogistic <- glm(formula = Forest ~ .,data = TreecoverBinarizedFormatted$GoodDat,family = "binomial")$coefficients


nSimsPerSetting <- 2000
SimsResultPath <- paste0(getwd(),"/Results/LowerNumLabels/")
#nVec <- c(50,100,150,200,250)
nVec <- c(250,200,150,100,50)
ExpsToRun <- c(2,3,4,6) #Change this to set which experiments are being run


#Experiment #1
if(1 %in% ExpsToRun){
  
  print("Experiment 1")
  
  
  SimsAlphaFoldLogReg <- list()
  for(nIdx in 1:length(nVec)){
    print(paste0("n= ",nVec[nIdx]))
    
    AlphaFoldSimsCurr <- runSimsRealData(FormattedDatInp = AlphaFoldFormatted,nSims = nSimsPerSetting,N_sampInp = 7500,nLabelTargInp = nVec[nIdx],regTypeInp = "logistic",outcomeVarName = "IDR")
    
    SimsAlphaFoldLogReg[[nIdx]] <- ExtractSummaryTableSingleExp(SimResultsInp = AlphaFoldSimsCurr$SimResultList,CoefNamesInp = CoefNameMatAlphaFold,nInp   = nVec[nIdx],TrueBeta = TrueBetaAlphaFoldLogistic, MethodNameMat = MethodNameMatPTD)
  }
  
  save(SimsAlphaFoldLogReg,nSimsPerSetting, file = paste0(SimsResultPath,"/AlphaFoldSimLessLabels.RData"))
}

#Experiment #2
if(2 %in% ExpsToRun){
  
  print("Experiment 2")
  
  
  SimsHousingPriceReg <- list()
  for(nIdx in 1:length(nVec)){
    print(paste0("n= ",nVec[nIdx]))
    HousingPriceSimsCurr <- runSimsRealData(FormattedDatInp = HousingPriceFormatted,nSims = nSimsPerSetting,N_sampInp = 5000,nLabelTargInp = nVec[nIdx],regTypeInp = "linear",outcomeVarName = "HousingPrice")
    
    SimsHousingPriceReg[[nIdx]] <- ExtractSummaryTableSingleExp(SimResultsInp = HousingPriceSimsCurr$SimResultList,CoefNamesInp = CoefNameMatHousingPrice,nInp   = nVec[nIdx],TrueBeta = TrueBetaHousingPriceLM, MethodNameMat = MethodNameMatPTD)
  }
  
  save(SimsHousingPriceReg,nSimsPerSetting, file = paste0(SimsResultPath,"/HousingPriceSimLessLabels.RData"))
}

#Experiment #3
if(3 %in% ExpsToRun){
  
  print("Experiment 3")
  
  SimsHousingPriceQuantileReg <- list()
  for(nIdx in 1:length(nVec)){
    print(paste0("n= ",nVec[nIdx]))
    HousingPriceSimsCurrQR <- runSimsRealData(FormattedDatInp = HousingPriceFormatted,nSims = nSimsPerSetting,N_sampInp = 5000,nLabelTargInp = nVec[nIdx],regTypeInp = "Quantile Regression",BootstrapSchemesInp = c("FullBootPercentile","QuickConvolution"),outcomeVarName = "HousingPrice",tauForQuantileReg = 0.5)
    
    
    
    SimsHousingPriceQuantileReg[[nIdx]] <- ExtractSummaryTableSingleExp(SimResultsInp = HousingPriceSimsCurrQR$SimResultList,CoefNamesInp = CoefNameMatHousingPrice,nInp   = nVec[nIdx],TrueBeta = TrueBetaHousingPriceQR, MethodNameMat = MethodNameMatPTD)
  }
  
  save(SimsHousingPriceQuantileReg,nSimsPerSetting, file = paste0(SimsResultPath,"/HousingPriceQuantileRegSimLessLabels.RData"))
}





if(4 %in% ExpsToRun){
  #Experiment #4
  print("Experiment 4")
  
  SimsTreecoverReg <- list()
  for(nIdx in 1:length(nVec)){
    print(paste0("n= ",nVec[nIdx]))
    TreecoverSimsCurr <- runSimsRealData(FormattedDatInp = TreecoverFormatted,nSims = nSimsPerSetting,N_sampInp = 5000,nLabelTargInp = nVec[nIdx],regTypeInp = "linear",outcomeVarName = "Treecover")
    
    SimsTreecoverReg[[nIdx]] <- ExtractSummaryTableSingleExp(SimResultsInp = TreecoverSimsCurr$SimResultList,CoefNamesInp = CoefNameMatTreecover,nInp   = nVec[nIdx],TrueBeta = TrueBetaTreecoverLM, MethodNameMat = MethodNameMatPTD)
  }

save(SimsTreecoverReg,nSimsPerSetting, file = paste0(SimsResultPath,"/TreecoverSimLessLabels.RData"))
}

#Experiment #5
if(5 %in% ExpsToRun){
  print("Experiment 5")
  
  SimsTreecoverClusteredReg <- list()
  for(nIdx in 1:length(nVec)){
    print(paste0("n= ",nVec[nIdx]))
    TreecoverSimsClustCurr <- runSimsRealData(FormattedDatInp = TreecoverFormatted,nSims = nSimsPerSetting,N_sampInp = 10000,nLabelTargInp = nVec[nIdx],regTypeInp = "linear",outcomeVarName = "Treecover",clusteredLabelling = T,useClusterIDs = T,BootstrapSchemesInp = c("QuickConvolution","FullBootPercentile"))
    SimsTreecoverClusteredReg[[nIdx]] <- ExtractSummaryTableSingleExp(SimResultsInp = TreecoverSimsClustCurr$SimResultList,CoefNamesInp = CoefNameMatTreecover,nInp   = nVec[nIdx],TrueBeta = TrueBetaTreecoverLM, MethodNameMat = MethodNameMatPTD)
  }
  
  save(SimsTreecoverClusteredReg,nSimsPerSetting, file = paste0(SimsResultPath,"/TreecoverClusteredSimLessLabels.RData"))
}



#Experiment #6
if(6 %in% ExpsToRun){
  print("Experiment 6")
  
  
  simsLowerLabelsForestLogReg <- list()
  for(nIdx in 1:length(nVec)){
    print(paste0("n= ",nVec[nIdx]))
    TreecoverBinarizedSimsCurr <- runSimsRealData(FormattedDatInp = TreecoverBinarizedFormatted,nSims = nSimsPerSetting,N_sampInp = 8000,nLabelTargInp = nVec[nIdx],regTypeInp = "logistic",BootstrapSchemesInp = c("FullBootPercentile","QuickConvolution"),outcomeVarName = "Forest")
    
    simsLowerLabelsForestLogReg[[nIdx]] <- ExtractSummaryTableSingleExp(SimResultsInp = TreecoverBinarizedSimsCurr$SimResultList,CoefNamesInp = CoefNameMatTreecover,nInp   = nVec[nIdx],TrueBeta = TrueBetaForestedLogistic, MethodNameMat = MethodNameMatPTD)
  }
  
  save(simsLowerLabelsForestLogReg,nSimsPerSetting, file = paste0(SimsResultPath,"/ForestLogRegSimLessLabels.RData"))

}
