#setwd("~/AdditionalInvestigations")




RandomSubsetFixedSize <- function(FormattedDatInp,N_samp){
  
  FormattedDatSubset <- list()
  
  
  MSuper <- nrow(FormattedDatInp$GoodDat)
  idxSubsample <- sample(1:MSuper,size = N_samp,replace = F) #If you sample with replacement simulations should have more precise coverage
  
  
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




library(dplyr)
SoftwareImplementationPath <- paste0(dirname(getwd()),"/PTD_Boot_Implementation/")
source(paste0(SoftwareImplementationPath,"PTDBootModularized.R"))
source(paste0(SoftwareImplementationPath,"BootstrapAndCalcEsts.R"))
source(paste0(SoftwareImplementationPath,"EstBetaGammaCalib_andVCOV.glm.R"))
source(paste0(SoftwareImplementationPath,"FormatAndFitReg.R"))


#implementing PPI++ (npMod and ppiMod are python modules. ppiMod is the ppi_py module)
library(reticulate)
npModule <- import("numpy")
#py_install("ppi-python",pip = TRUE) #Install ppi_py
py_require(c("ppi-python"))
ppiModule <- import("ppi_py")

RunPPIPlusPlusglm <- function(GoodDat,ProxyDat,RegType,OutcomeVarName,alphaInp2,npMod,ppiMod){
  #Formatting data
  idxComp <- complete.cases(GoodDat)
  YAll <- GoodDat[[OutcomeVarName]]
  YhatAll <- ProxyDat[[OutcomeVarName]]
  XAllNoInt <- ProxyDat #There should be no Error-in-covariates
  XAllNoInt[[OutcomeVarName]] <- NULL #dropping outcome variable
  XAll <- cbind(rep(1,nrow(XAllNoInt)),XAllNoInt) #Adding intercept
  CovariateNames <- c("(Intercept)",names(XAllNoInt))
  
  
  if(RegType %in% c("linear","Linear","Gaussian","gaussian","Normal","normal","OLS","ols")){
    ptEstFunct <- ppiMod$ppi_ols_pointestimate
    CIFunct <- ppiMod$ppi_ols_ci
  } else if(RegType %in% c("logistic","Logistic","logit")){
    ptEstFunct <- ppiMod$ppi_logistic_pointestimate
    CIFunct <- ppiMod$ppi_logistic_ci
  }
  
  ptEst <- ptEstFunct(X = npMod$array(as.matrix(XAll[idxComp, ])),
                      Y           = npMod$array(as.numeric(YAll[idxComp])),
                      Yhat        = npMod$array(as.numeric(YhatAll[idxComp])),
                      X_unlabeled = npMod$array(as.matrix(XAll[!idxComp, ])),
                      Yhat_unlabeled = npMod$array(as.numeric(YhatAll[!idxComp]))
  )
  
  CIcurr <- CIFunct(X = npMod$array(as.matrix(XAll[idxComp, ])),
                    Y           = npMod$array(as.numeric(YAll[idxComp])),
                    Yhat        = npMod$array(as.numeric(YhatAll[idxComp])),
                    X_unlabeled = npMod$array(as.matrix(XAll[!idxComp, ])),
                    Yhat_unlabeled = npMod$array(as.numeric(YhatAll[!idxComp])),
                    alpha       = alphaInp2
  )
  
  ResultMatOut <- cbind(ptEst,CIcurr[[1]],CIcurr[[2]])
  rownames(ResultMatOut) <- CovariateNames
  colnames(ResultMatOut) <- c("Estimate","CI_lb","CI_ub")
  return(ResultMatOut)
}



library(doParallel)
runSimsRealDataWithPPI <- function(FormattedDatInp2,N_sampInp,nLabelTargInp,nSims=100,regTypeInp="linear",outcomeVarName="Y",BootstrapSchemesInp=c("FullBootPercentile","QuickConvolution","CLTBased"),tuningSchemesInp=c("Optimal","Diagonal","None"),alphaInp=0.1,tauForQuantileReg=NULL,runPPI=T,npModInp=NULL,ppiModInp=NULL,ncores=78){
  
  registerDoParallel(ncores)
  
  
  currTime <- proc.time()[3]
  startTime <- currTime
  
  writeLines("", "progress.log")
  list()
  currTime <- proc.time()[3]
  startTime <- currTime
  SimResultList <- foreach(i = 1:nSims) %dopar% {
    
    cat(paste0("Starting Simulation ", i, " out of ", nSims,
               " | ", format(Sys.time(), "%H:%M:%S"), "\n"),
        file = "progress.log", append = TRUE)
    
    
    
    FormattedDatSubsamp <- RandomSubsetFixedSize(FormattedDatInp = FormattedDatInp2,N_samp = N_sampInp)
    
    CurrLabelling <- LabelRandomSubset(GoodDatInp = FormattedDatSubsamp$GoodDat,
                                       ProxyDatInp = FormattedDatSubsamp$ProxyDat,
                                       piLabelInpUsnc = FormattedDatSubsamp$PiLabelUnsc,
                                       nLabelTarget = nLabelTargInp)
    
    
    currSimResultsAllMethods <- list()
    suppressWarnings(currSimResultsAllMethods$PTDSuite <- PTDBootModularized(ProxyDat = FormattedDatSubsamp$ProxyDat,GoodDat = CurrLabelling$GoodDatWMissing,PiLabel = CurrLabelling$piLabelRescaled ,BootstrapScheme=BootstrapSchemesInp,nBootInference = 2000,nBootTune = 0,TuningScheme = tuningSchemesInp,RegType = regTypeInp,OutcomeVarName = outcomeVarName,alpha = alphaInp,tauQuantReg = tauForQuantileReg)$CIList) 
    
    ComparisonMethodsCurr <- list()
    if(runPPI){
      ComparisonMethodsCurr$PPIplusPlus <- RunPPIPlusPlusglm(GoodDat = CurrLabelling$GoodDatWMissing,ProxyDat = FormattedDatSubsamp$ProxyDat,RegType = regTypeInp,OutcomeVarName = outcomeVarName,alphaInp2 = alphaInp,npMod = npModInp ,ppiMod = ppiModInp)
      
    }
    
    currSimResultsAllMethods$ComparisonMethodsSuite <- ComparisonMethodsCurr
    
    
    currSimResultsAllMethods 
  }
  
  file.remove("progress.log")
  stopImplicitCluster()
  
  
  
  SimResultListPTDAndClass <- list()
  ComparisonMethods <- list()
  for(i in 1:nSims){
    SimResultListPTDAndClass[[i]] <- SimResultList[[i]]$PTDSuite
    ComparisonMethods[[i]] <- SimResultList[[i]]$ComparisonMethodsSuite
  }
  
  totalTime <- (proc.time()[3]-startTime)/60
  return(list(SimResultListPTDAndClass=SimResultListPTDAndClass,ComparisonMethods=ComparisonMethods,N_sampInp=N_sampInp,nLabelTargInp=nLabelTargInp,RegTypeUsed=regTypeInp,alphaUsed=alphaInp,tauForQuantileRegUsed=tauForQuantileReg,totalTime=totalTime)) 
}



ExtractSummaryTableSingleExp <- function(SimResultsInp,CoefNamesInp,MethodNameMat,TrueBeta,ErrorParam){
  
  MethodNameVec <-  MethodNameMat[,1] 
  MethodNameFormatted <- MethodNameMat[,2] 
  df1TuneMethodList <- list()
  for(methId in 1:length(MethodNameVec)){
    nCoeff <- nrow(CoefNamesInp)
    
    if(MethodNameVec[methId] %in% c("PPIplusPlus","MultipleImputation")){
      ResultListCurrMethod <- SimResultsInp$ComparisonMethods
    } else {
      ResultListCurrMethod <- SimResultsInp$SimResultListPTDAndClass
    }
    
    CIWidthsTabCurr <- matrix(NA,ncol = length(ResultListCurrMethod),nrow = nCoeff)
    CoversMat <- matrix(NA,ncol = length(ResultListCurrMethod),nrow = nCoeff)
    
    for(j in 1:length(ResultListCurrMethod)){
      
      CIsCurrSim <- ResultListCurrMethod[[j]][[MethodNameVec[methId]]]
      
      
      if((sum(rownames(CIsCurrSim) != CoefNamesInp[,1])==0) | (sum(rownames(CIsCurrSim)  != names(TrueBeta))==0) ){
        CIWidthsTabCurr[,j] <- CIsCurrSim[,'CI_ub']-CIsCurrSim[,'CI_lb']
        CoversMat[,j] <- I((TrueBeta<CIsCurrSim[,'CI_ub']) & (TrueBeta > CIsCurrSim[,'CI_lb']))
      } else{print("Warning:Coefficient name misalignment")}
    }
    
    CIWidths_mean <- rowMeans(CIWidthsTabCurr)
    CIWidths_sd <- rep(NA,nCoeff)
    for(i in 1:nCoeff){CIWidths_sd[i] <- sd(CIWidthsTabCurr[i,])}
    EmpCoverage <- rowMeans(CoversMat)
    
    dfOutCurr <- data.frame(CoefNamesInp[,2],rep(MethodNameFormatted[methId],nCoeff),
                            rep(ErrorParam,nCoeff),
                            CIWidths_mean,CIWidths_sd,EmpCoverage,TrueBeta)
    names(dfOutCurr) <- c("Coefficient","Method","ErrorParam","CIWidths_mean","CIWidths_sd","Coverage","TrueBeta")
    df1TuneMethodList[[methId]] <- dfOutCurr
  }
  dfOut <- df1TuneMethodList[[1]]
  for(methId in 2:length(MethodNameVec)){dfOut <- rbind(dfOut,df1TuneMethodList[[methId]])}
  return(dfOut)
}



#Running 3 Experiments

#1) With Error-in-Outcome Housing Price Linear regression
#2) With Error-in-Outcome Treecover Linear regression
#3) With Error-in-Outcome Forest Logistic regression

#We will standardize all the datasets so that trace tuning doesn't priortize any particular coefficient.

#In each setting we translate the proxies and observe whether PPI++ and PTD degrades.


#Format datasets to be Error-in-outcome only (no missing covariates) and standardized

#loading datasets
load(paste0(dirname(getwd()),"/Datasets/","HousingPriceFormattedWithProxyHP.RData"))
load(paste0(dirname(getwd()),"/Datasets/","TreecoverBinarizedFormatted.RData"))
load(paste0(dirname(getwd()),"/Datasets/","TreecoverFormatted.RData"))



#Function to remove error in covariates
RemoveErrorInCovs <- function(dfInp,OutcomeVarName){
  dfOut <- dfInp
  dfOut$ProxyDat <- dfInp$GoodDat
  dfOut$ProxyDat[[OutcomeVarName]] <- dfInp$ProxyDat[[OutcomeVarName]]
  return(dfOut)
}

#Function to standardize dataset according to mean and sd of good data
stdzDataset <- function(dfInp2){
  dfOut2 <- dfInp2
  varsDf <- names(dfInp2$GoodDat)
  for(j in 1:length(varsDf)){
    currVar <- varsDf[j]
    currVarVecGood <- dfInp2$GoodDat[[currVar]]
    if(!is.logical(currVarVecGood)){ #Don't standardize binary variables
meanCurr <- mean(currVarVecGood)
sdCurr <- sd(currVarVecGood)
dfOut2$GoodDat[[currVar]] <- (currVarVecGood-meanCurr)/sdCurr
dfOut2$ProxyDat[[currVar]] <- (dfInp2$ProxyDat[[currVar]]-meanCurr)/sdCurr
} else{
  print(paste0(currVar," is a logical. Not standardizing it."))
}
}
return(dfOut2)
}


HousingPriceErrorInOutcomeStdz <- stdzDataset(RemoveErrorInCovs(HousingPriceFormattedWithProxyHP,OutcomeVarName="HousingPrice"))

TreecoverErrorInOutcomeStdz <- stdzDataset(RemoveErrorInCovs(TreecoverFormatted,OutcomeVarName="Treecover"))

ForestBinarizedErrorInOutcomeStdz <- stdzDataset(RemoveErrorInCovs(TreecoverBinarizedFormatted,OutcomeVarName="Forest"))

#Checking that only the outcome variable has prediction errors
colMeans(HousingPriceErrorInOutcomeStdz$ProxyDat==HousingPriceErrorInOutcomeStdz$GoodDat)
colMeans(TreecoverErrorInOutcomeStdz$ProxyDat==TreecoverErrorInOutcomeStdz$GoodDat)
colMeans(ForestBinarizedErrorInOutcomeStdz$ProxyDat==ForestBinarizedErrorInOutcomeStdz$GoodDat)

#Removing original dataframes
rm(TreecoverFormatted,TreecoverBinarizedFormatted,HousingPriceFormattedWithProxyHP)





#Function that shifts the mean of the predictions (either through translation in continuous outcome case or flipping a random fraction of the predictions)

ShiftPredictionMean <- function(dfInp,OutcomeVarName,tInp){
  OutcomeVarProxyOriginal <- dfInp$ProxyDat[[OutcomeVarName]]
  dfOut <- dfInp
  if(is.logical(OutcomeVarProxyOriginal)){
    Ntot <- length(OutcomeVarProxyOriginal)
    OutcomeVarProxyNew <- OutcomeVarProxyOriginal
    if(tInp>0){ #If outcome variable is binary and want to increase mean
      idxProxy0 <- which(OutcomeVarProxyOriginal==FALSE)
      idxTurnTo1 <- sample(idxProxy0,size = round(Ntot*tInp) ,replace = F)
      OutcomeVarProxyNew[idxTurnTo1] <- TRUE # flip appropriate # of predictions from 0 to 1
    } else if(tInp<0){ #If outcome variable is binary and want to decrease mean
      idxProxy1 <- which(OutcomeVarProxyOriginal==TRUE)
      idxProxy0 <- sample(idxProxy1,size = round(Ntot*abs(tInp)) ,replace = F)
      OutcomeVarProxyNew[idxProxy0] <- FALSE # flip appropriate # of predictions from 1 to 0
    }
    dfOut$ProxyDat[[OutcomeVarName]] <- OutcomeVarProxyNew 
  } else{ #If outcome variable is continuous just translate the predictions directly
    dfOut$ProxyDat[[OutcomeVarName]] <- tInp+OutcomeVarProxyOriginal 
  }
  return(dfOut)
}



#Specifying input parameters
MethodNameMatWithPPI <- cbind(c("QuickConvolution_Diagonal","QuickConvolution_ScalarTrace","PPIplusPlus","classical"),c("PTD (Optimal Diagonal Tuning)","PTD (Scalar Tuning)","PPI++","Classical"))

CoefNameMatHousingPrice <- cbind(c("(Intercept)","Income","Nightlights","RoadLength"),
                                 c("Intercept","Income","Nightlights","RoadLength"))

CoefNameMatTreecover <- cbind(c("(Intercept)","Elevation","Population"),
                              c("Intercept","Elevation","Population"))

#Calculating True Regression Coefficients
TrueBetaHousingPriceLM <- lm(formula = HousingPrice ~ .,data = HousingPriceErrorInOutcomeStdz$GoodDat)$coefficients

TrueBetaTreecoverLM <- lm(formula = Treecover ~ .,data = TreecoverErrorInOutcomeStdz$GoodDat)$coefficients

TrueBetaForestedLogistic <- glm(formula = Forest ~ .,data = ForestBinarizedErrorInOutcomeStdz$GoodDat,family = "binomial")$coefficients

nSimsPerSetting <- 2000
SimsResultPath <- paste0(getwd(),"/Results/PPI++Comp/")



tVec <- seq(-3,3,by = 0.5)

SimsHousingPrice <- list()

for(t in 1:length(tVec)){
  print(tVec[t])
  FormattedDatModif <- ShiftPredictionMean(dfInp = HousingPriceErrorInOutcomeStdz,tInp = tVec[t],OutcomeVarName = "HousingPrice")
  
  SimsHousingPriceCurr <- runSimsRealDataWithPPI(FormattedDatInp2 = FormattedDatModif,N_sampInp = 5000,nLabelTargInp = 500,nSims = nSimsPerSetting,regTypeInp = "linear",outcomeVarName = "HousingPrice",BootstrapSchemesInp = "QuickConvolution",tuningSchemesInp = c("Diagonal","ScalarTrace"),alphaInp = 0.1,runPPI = TRUE,npModInp = npModule,ppiModInp = ppiModule)
  
  SimsHousingPrice[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = SimsHousingPriceCurr,CoefNamesInp = CoefNameMatHousingPrice,TrueBeta = TrueBetaHousingPriceLM ,MethodNameMat = MethodNameMatWithPPI,ErrorParam  = tVec[t])
}


save(SimsHousingPrice,nSimsPerSetting, file = paste0(SimsResultPath,"/IncreaseTranslationErrorPPI++HousingPriceQuickBoot.RData"))





tVec <- seq(-3,3,by = 0.5)

SimsTreecover <- list()

for(t in 1:length(tVec)){
  print(tVec[t])
  FormattedDatModif <- ShiftPredictionMean(dfInp = TreecoverErrorInOutcomeStdz,tInp = tVec[t],OutcomeVarName = "Treecover")
  
  SimsTreecoverCurr <- runSimsRealDataWithPPI(FormattedDatInp2 = FormattedDatModif,N_sampInp = 5000,nLabelTargInp = 500,nSims = nSimsPerSetting,regTypeInp = "linear",outcomeVarName = "Treecover",BootstrapSchemesInp = "QuickConvolution",tuningSchemesInp = c("Diagonal","ScalarTrace"),alphaInp = 0.1,runPPI = TRUE,npModInp = npModule,ppiModInp = ppiModule)
  
  SimsTreecover[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = SimsTreecoverCurr,CoefNamesInp = CoefNameMatTreecover,TrueBeta = TrueBetaTreecoverLM,MethodNameMat = MethodNameMatWithPPI,ErrorParam  = tVec[t])
}


save(SimsTreecover,nSimsPerSetting, file = paste0(SimsResultPath,"/IncreaseTranslationErrorPPI++TreecoverQuickBoot.RData"))




tVec <- seq(-3,3,by = 0.5)/10

SimsForestCoverLogReg <- list()

for(t in 1:length(tVec)){
  print(tVec[t])
  FormattedDatModif <- ShiftPredictionMean(dfInp = ForestBinarizedErrorInOutcomeStdz,tInp = tVec[t],OutcomeVarName = "Forest")
  
  SimsForestCurr <- runSimsRealDataWithPPI(FormattedDatInp2 = FormattedDatModif,N_sampInp = 8000,nLabelTargInp = 1000,nSims = nSimsPerSetting,regTypeInp = "logistic",outcomeVarName = "Forest",BootstrapSchemesInp = "QuickConvolution",tuningSchemesInp = c("Diagonal","ScalarTrace"),alphaInp = 0.1,runPPI = TRUE,npModInp = npModule,ppiModInp = ppiModule)
  
  SimsForestCoverLogReg[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = SimsForestCurr,CoefNamesInp = CoefNameMatTreecover,TrueBeta = TrueBetaForestedLogistic,MethodNameMat = MethodNameMatWithPPI,ErrorParam  = tVec[t])
}


save(SimsForestCoverLogReg,nSimsPerSetting, file = paste0(SimsResultPath,"/IncreaseTranslationErrorPPI++ForestLogRegQuickBoot.RData"))
