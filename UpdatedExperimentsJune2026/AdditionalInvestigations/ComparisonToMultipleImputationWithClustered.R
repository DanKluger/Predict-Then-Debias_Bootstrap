#setwd("~/AdditionalInvestigations")
library(quantreg)



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

#Implement multiple imputation
library(sandwich)
library(mice)
library(mitools)

RunMultipleImp <- function(GoodDat,ProxyDat,VarsToImpute,alphaInp2=0.1,mImpute=100,RegSpecId,ClusterID=NULL){
  
  #Formatting
  
  VarsWithMissing <- GoodDat[,which(names(GoodDat) %in% VarsToImpute)]
  VarNames <- names(ProxyDat)
  ProxyDatVarNames <- VarNames
  for(j in 1:length(VarNames)){
    if(VarNames[j] %in% VarsToImpute){
      ProxyDatVarNames[j] <- paste0("proxy_",VarNames[j])
    }
  }
  names(ProxyDat) <- ProxyDatVarNames
  
  dfMI <- cbind(VarsWithMissing,ProxyDat)
 # dfMI$ClusterID <- ClusterID
  if(length(VarsToImpute)==1){names(dfMI)[1] <- VarsToImpute}
  
  # Set up the imputation method: only Y gets imputed
  meth <- make.method(dfMI)
  meth[] <- ""          # don't impute anything by default
  for(j in 1:length(VarsToImpute)){
    curVar <- VarsToImpute[j]
    if(is.logical(dfMI[[curVar]])){
      meth[curVar] <- "logreg"   # Bayesian logistic regression for binary variables
    } else if(is.numeric(dfMI[[curVar]])){
      meth[curVar] <- "norm"   # Bayesian linear regression for continuous variables
    }
  }
  
  # Set up the predictor matrix: which variables predict missing values in the imputation model
  pred <- make.predictorMatrix(dfMI)
  pred[,] <- 0                          # clear everything
  pred[VarsToImpute, ProxyDatVarNames] <- 1 
  
  
  # Run MICE
  imp <- mice(dfMI, method = meth, predictorMatrix = pred, m = mImpute, maxit = 30, print = FALSE)
  
  # Fit the regression on each imputed dataset
  if(RegSpecId=="HousingPriceReg"){fit <- with(imp, lm(HousingPrice ~ Income + Nightlights + RoadLength))}
  if(RegSpecId=="HousingPriceQuantReg"){fit <- with(imp, rq(HousingPrice ~ Income + Nightlights + RoadLength,tau=0.5))}
  if(RegSpecId=="TreecoverReg"){fit <- with(imp, lm(Treecover ~ Elevation + Population))}
  if(RegSpecId=="ForestedLogReg"){fit <- with(imp, glm(Forest ~ Elevation + Population,family = "binomial"))}

    
  if(RegSpecId!="TreecoverClusteredReg"){
    # Pool the results using Rubin's rules
    pooled <- pool(fit)
    ResultsMatFull <- summary(pooled, conf.int = TRUE, conf.level = 1-alphaInp2)
  
    #Extract and return results of interest
  
    # Find columns that match the "number %" pattern
    pct_cols <- grep("^[0-9.]+ %$", names(ResultsMatFull), value = TRUE)
  
    # Extract the numeric values to find the lowest
    pct_values <- as.numeric(sub(" %", "", pct_cols))
  
    # Get the column name with the lowest and highest percentage
    low_col <- pct_cols[which.min(pct_values)]
    high_col <- pct_cols[which.max(pct_values)]
  
    ResultsMatOut <- as.matrix(cbind(ResultsMatFull$estimate,ResultsMatFull[[low_col]],ResultsMatFull[[high_col]]))
    rownames(ResultsMatOut) <- ResultsMatFull$term
  } else if(RegSpecId=="TreecoverClusteredReg"){
    # Manual fit + clustered vcov on each imputed dataset
    betas <- vector("list", imp$m)
    vcovs <- vector("list", imp$m)
    for (i in seq_len(imp$m)) {
      d_i        <- complete(imp, i)
      mod_i      <- lm(Treecover ~ Elevation + Population, data = d_i)
      betas[[i]] <- coef(mod_i)
      vcovs[[i]] <- vcovCL(mod_i, cluster = ClusterID,sandwich = T)
    }
    
    # Pool using Rubin's rules via mitools
    pooled_mi <- MIcombine(results = betas, variances = vcovs)
    
    est   <- coef(pooled_mi)
    se    <- sqrt(diag(vcov(pooled_mi)))
    df    <- pooled_mi$df
    tcrit <- qt(1 - alphaInp2 / 2, df = df)
    
    ResultsMatOut <- as.matrix(cbind(est,est - tcrit * se,est + tcrit * se))
    rownames(ResultsMatOut) <- names(est)
  }
  
  colnames(ResultsMatOut) <- c("Estimate","CI_lb","CI_ub")
  return(ResultsMatOut)
}

library(doParallel)
runSimsRealDataWithMI <- function(FormattedDatInp2,N_sampInp,nLabelTargInp,nSims=100,regTypeInp="linear",outcomeVarName="Y",BootstrapSchemesInp=c("FullBootPercentile","QuickConvolution","CLTBased"),tuningSchemesInp=c("Optimal","Diagonal","None"),alphaInp=0.1,tauForQuantileReg=NULL,runMI=T,RegSpecIdInp=NULL,ncores=70,clusteredLabelling=F,useClusterIDs=F){
  
  
  SimResultListPTDAndClass <- list()
  ComparisonMethods <- list()
  varsSometimesMissing <- vector()
  varNamesAll <- names(FormattedDatInp2$GoodDat)
  for(j in 1:length(varNamesAll)){
    ident <- identical(FormattedDatInp2$GoodDat[[varNamesAll[j]]],FormattedDatInp2$ProxyDat[[varNamesAll[j]]])
    if(!ident){ varsSometimesMissing <- c(varsSometimesMissing,varNamesAll[j])}
  }
  print(varsSometimesMissing)
  
  
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
    
    
    
    currSimResultsAllMethods <- list()
    
    suppressWarnings(currSimResultsAllMethods$PTDSuite <- PTDBootModularized(ProxyDat = FormattedDatSubsamp$ProxyDat,GoodDat = CurrLabelling$GoodDatWMissing,PiLabel = CurrLabelling$piLabelRescaled,clusterID=ClustersUse,BootstrapScheme=BootstrapSchemesInp,nBootTune = 0,TuningScheme = tuningSchemesInp,RegType = regTypeInp,OutcomeVarName = outcomeVarName,alpha = alphaInp,tauQuantReg = tauForQuantileReg)$CIList) 
    
    ComparisonMethodsCurr <- list()
    
    if(runMI){
      ComparisonMethodsCurr$MultipleImputation <- RunMultipleImp(GoodDat = CurrLabelling$GoodDatWMissing,ProxyDat = FormattedDatSubsamp$ProxyDat,ClusterID=ClustersUse,VarsToImpute = varsSometimesMissing,alphaInp2 = alphaInp,RegSpecId = RegSpecIdInp) 
    }
    
    currSimResultsAllMethods$ComparisonMethodsSuite <- ComparisonMethodsCurr
    currSimResultsAllMethods
  }
  
  
  file.remove("progress.log")
  stopImplicitCluster()
  
  totalTime <- (proc.time()[3]-startTime)/60
  
  

  for(i in 1:nSims){
    SimResultListPTDAndClass[[i]] <- SimResultList[[i]]$PTDSuite
    ComparisonMethods[[i]] <- SimResultList[[i]]$ComparisonMethodsSuite
  }
  
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
    
    dfOutCurr <- data.frame(CoefNamesInp[,2],rep(MethodNameFormatted[methId],nCoeff),rep(ErrorParam,nCoeff),
                            CIWidths_mean,CIWidths_sd,EmpCoverage,TrueBeta)
    names(dfOutCurr) <- c("Coefficient","Method","ErrorParam","CIWidths_mean","CIWidths_sd","Coverage","TrueBeta")
    df1TuneMethodList[[methId]] <- dfOutCurr
  }
  dfOut <- df1TuneMethodList[[1]]
  for(methId in 2:length(MethodNameVec)){dfOut <- rbind(dfOut,df1TuneMethodList[[methId]])}
  return(dfOut)
}



#Function that raises a predictor variable to a certain power

ExponentiateProxy <- function(dfInp,varToExp,lambda){
  ProxyOriginalCurrVar <- dfInp$ProxyDat[[varToExp]]
  ProxyExpUnnormalized <- ProxyOriginalCurrVar^lambda
  dfOut <- dfInp
  dfOut$ProxyDat[[varToExp]] <- ProxyExpUnnormalized*mean(ProxyOriginalCurrVar)/mean(ProxyExpUnnormalized)
  return(dfOut)
}

#Load datasets
load(paste0(dirname(getwd()),"/Datasets/","HousingPriceFormatted.RData"))
load(paste0(dirname(getwd()),"/Datasets/","TreecoverBinarizedFormatted.RData"))
load(paste0(dirname(getwd()),"/Datasets/","TreecoverFormatted.RData"))



MethodNameMatWithMI <- cbind(c("QuickConvolution_Diagonal","MultipleImputation","classical"),c("PTD (Optimal Diagonal Tuning)","Multiple Imputation","Classical"))

CoefNameMatTreecover <- cbind(c("(Intercept)","Elevation","Population"),
                              c("Intercept","Elevation","Population"))


CoefNameMatHousingPrice <- cbind(c("(Intercept)","Income","Nightlights","RoadLength"),
                                 c("Intercept","Income","Nightlights","RoadLength"))



#Calculating True Regression Coefficients
TrueBetaHousingPriceLM <- lm(formula = HousingPrice ~ .,data = HousingPriceFormatted$GoodDat)$coefficients
TrueBetaHousingPriceQuantReg <- rq(formula = HousingPrice ~ .,data = HousingPriceFormatted$GoodDat,tau=0.5)$coefficients


TrueBetaTreecoverLM <- lm(formula = Treecover ~ .,data = TreecoverFormatted$GoodDat)$coefficients

TrueBetaForestedLogistic <- glm(formula = Forest ~ .,data = TreecoverBinarizedFormatted$GoodDat,family = "binomial")$coefficients


nSimsPerSetting <- 2000
lambdaVec <- seq(1,5) 
SimsResultPath <- paste0(getwd(),"/Results/MultipleImputationComp/")
ExpsToRun <- c(2,3,4,5,6)


library(dplyr)

if(2 %in% ExpsToRun){
  SimsHousingPrice <- list()
  
  for(t in 1:length(lambdaVec)){
    print(paste0("lambda=",lambdaVec[t]))
    FormattedDatModif <- ExponentiateProxy(dfInp = HousingPriceFormatted,lambda = lambdaVec[t],varToExp = "RoadLength")
    
    SimsHousingPriceCurr <- runSimsRealDataWithMI(FormattedDatInp2 = FormattedDatModif,N_sampInp = 5000,nLabelTargInp = 500,nSims = nSimsPerSetting,regTypeInp = "linear",outcomeVarName = "HousingPrice",BootstrapSchemesInp = "QuickConvolution",tuningSchemesInp = "Diagonal",alphaInp = 0.1,runMI = TRUE,RegSpecIdInp = "HousingPriceReg")
    
    SimsHousingPrice[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = SimsHousingPriceCurr,CoefNamesInp = CoefNameMatHousingPrice,TrueBeta = TrueBetaHousingPriceLM ,MethodNameMat = MethodNameMatWithMI,ErrorParam  = lambdaVec[t])
  }
  
  
  save(SimsHousingPrice,nSimsPerSetting, file = paste0(SimsResultPath,"/ExponentiateProxyHousingPriceWithMI.RData"))
}

if(3 %in% ExpsToRun){
  SimsHousingPriceQuantReg <- list()
  
  for(t in 1:length(lambdaVec)){
    print(paste0("lambda=",lambdaVec[t]))
    FormattedDatModif <- ExponentiateProxy(dfInp = HousingPriceFormatted,lambda = lambdaVec[t],varToExp = "RoadLength")
    
    SimsHousingPriceCurr <- runSimsRealDataWithMI(FormattedDatInp2 = FormattedDatModif,N_sampInp = 5000,nLabelTargInp = 1000,nSims = nSimsPerSetting,regTypeInp = "Quantile",outcomeVarName = "HousingPrice",BootstrapSchemesInp = "QuickConvolution",tuningSchemesInp = "Diagonal",alphaInp = 0.1,runMI = TRUE,RegSpecIdInp = "HousingPriceQuantReg",tauForQuantileReg = 0.5)
    
    SimsHousingPriceQuantReg[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = SimsHousingPriceCurr,CoefNamesInp = CoefNameMatHousingPrice,TrueBeta = TrueBetaHousingPriceQuantReg ,MethodNameMat = MethodNameMatWithMI,ErrorParam  = lambdaVec[t])
  }
  
  
  save(SimsHousingPriceQuantReg,nSimsPerSetting, file = paste0(SimsResultPath,"/ExponentiateProxyHousingPriceQuantRegWithMI.RData"))
}


if(4 %in% ExpsToRun){
  SimsTreecover <- list()
  
  for(t in 1:length(lambdaVec)){
    print(paste0("lambda=",lambdaVec[t]))
    FormattedDatModif <- ExponentiateProxy(dfInp = TreecoverFormatted,lambda = lambdaVec[t],varToExp = "Population")
    
    SimsTreecoverCurr <-  runSimsRealDataWithMI(FormattedDatInp2 = FormattedDatModif,N_sampInp = 5000,nLabelTargInp = 500,nSims = nSimsPerSetting,regTypeInp = "linear",outcomeVarName = "Treecover",BootstrapSchemesInp = "QuickConvolution",tuningSchemesInp = "Diagonal",alphaInp = 0.1,runMI = TRUE,RegSpecIdInp = "TreecoverReg")
    
    
    SimsTreecover[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = SimsTreecoverCurr,CoefNamesInp = CoefNameMatTreecover,TrueBeta = TrueBetaTreecoverLM,MethodNameMat = MethodNameMatWithMI,ErrorParam  = lambdaVec[t])
  }
  
  
  save(SimsTreecover,nSimsPerSetting, file = paste0(SimsResultPath,"/ExponentiateProxyTreecoverWithMI.RData"))
}


if(5 %in% ExpsToRun){
  SimsTreecoverClustered <- list()
  
  for(t in 1:length(lambdaVec)){
    print(paste0("lambda=",lambdaVec[t]))
    FormattedDatModif <- ExponentiateProxy(dfInp = TreecoverFormatted,lambda = lambdaVec[t],varToExp = "Population")
    
    SimsTreecoverCurr <-  runSimsRealDataWithMI(FormattedDatInp2 = FormattedDatModif,N_sampInp = 10000,nLabelTargInp = 1000,nSims = nSimsPerSetting,regTypeInp = "linear",outcomeVarName = "Treecover",BootstrapSchemesInp = "QuickConvolution",clusteredLabelling = T,useClusterIDs = T,tuningSchemesInp = "Diagonal",alphaInp = 0.1,runMI = TRUE,RegSpecIdInp = "TreecoverClusteredReg")
    SimsTreecoverClustered[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = SimsTreecoverCurr,CoefNamesInp = CoefNameMatTreecover,TrueBeta = TrueBetaTreecoverLM,MethodNameMat = MethodNameMatWithMI,ErrorParam  = lambdaVec[t])
  }
  
  
  save(SimsTreecoverClustered,nSimsPerSetting, file = paste0(SimsResultPath,"/ExponentiateProxyTreecoverClusteredWithMI.RData"))
}


if(6 %in% ExpsToRun){
  SimsForestCoverLogReg <- list()
  
  for(t in 1:length(lambdaVec)){
    print(paste0("lambda=",lambdaVec[t]))
    
    FormattedDatModif <- ExponentiateProxy(dfInp = TreecoverBinarizedFormatted,lambda = lambdaVec[t],varToExp = "Population")
    
    SimsForestCurr <- runSimsRealDataWithMI(FormattedDatInp2 = FormattedDatModif,N_sampInp = 8000,nLabelTargInp = 1000,nSims = nSimsPerSetting,regTypeInp = "logistic",outcomeVarName = "Forest",BootstrapSchemesInp = "QuickConvolution",tuningSchemesInp = "Diagonal",alphaInp = 0.1,runMI = TRUE,RegSpecIdInp = "ForestedLogReg")
    
    SimsForestCoverLogReg[[t]] <- ExtractSummaryTableSingleExp(SimResultsInp = SimsForestCurr,CoefNamesInp = CoefNameMatTreecover,TrueBeta = TrueBetaForestedLogistic,MethodNameMat = MethodNameMatWithMI,ErrorParam  = lambdaVec[t])
  }
  
  
  save(SimsForestCoverLogReg,nSimsPerSetting, file = paste0(SimsResultPath,"/ExponentiateProxyForestLogRegWithMI.RData"))
}

