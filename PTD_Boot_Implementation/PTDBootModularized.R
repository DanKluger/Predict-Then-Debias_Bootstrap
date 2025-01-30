##########################################################
# INPUTS
#  ProxyDat := An m x (p+1) data frame where all variables that are not widely available are replaced by their widely available proxies
#  GoodDat  := An m x (p+1) data frame with the actual data. Rows corresponding to data points that were not labelled should be set to NA
#  NOTE: ProxyDat and GoodDat should have the same column names 

#  PiLabel := A vector of length m giving the probability that a label was collected for each row in the proxy data.
#  NOTE: These probabilities need to be properly scaled such that they are probabilities (if they are merely proportional to the probabilities this can lead to issues)
#  If the labelled data is from the same distribution as the unlabelled data to insert a constant vector for PiLabel but it should be appropriately scaled

#  OutcomeVarName   := The name of the outcome variable in your data frame. This function will estimate the regression of that variable on the remaining variables.

#  RegType   := Use "linear" for linear regression, "logistic" for logistic regression,"Poisson" for Poisson regression, "Quantile" for Quantile regression (set tauQuantReg to be the quantile )

#  tauQuantReg := is the quantile of interest for quantile regression (only use this if you are running quantile regression)

#  nBootInference := The number of bootstrap draws B (the default is 2000. Larger B cannot hurt, but can lead to slower runtime)

#  alpha := Statistical significance level. The default is 0.05, but note that the paper used 0.1

#  BootstrapScheme :=   
#  Note that you can select multiple options in one run of this function
#  use  "CLTBased" for Algorithm 4 from the paper. This computes a confidence interval based on a plug in estimate for the asymptotic variance and the CLT (CLTBased CIs are currently available for GLM only and are not currently available for cluster sampling settings)
#  use "QuickConvolution" for Algorithm 3 from the paper. This an approach that takes a convolution of the normal approximation for the noncalibration term and a draw from the bootstrap distribution for the bias correction term
#  use "FullBootStudentized" for the slow (potentially 2nd order accurate) approach (studentized bootstrap is available for GLM only and not available for cluster sampling settings)
#  use "FullBasicBoot" for the slow basic bootstrap approach
#  use "FullBootPercentile" for Algorithm 2 from the paper. This is the percentile bootstrap approach

#  TuningScheme := character vector of tuning schemes to use. Note that you can select multiple options in one run of this function (at no major cost to the runtime)
#  "Optimal"  Estimates optimal tuning matrix
#  "Diagonal" Estimates optimal tuning matrix among diagonal tuning matrices
#  "None" sets Omega to the identity
#  "Zero" sets tuning matrix to zero (used for debugging purposes to see how well bootstrapping classical estimator works)

##########################################################



library(sandwich)
library(lmtest)

PTDBootModularized <-function(ProxyDat,GoodDat,PiLabel=NULL,clusterID=NULL,OutcomeVarName="Y",RegType="linear",BootstrapScheme="QuickConvolution",nBootTune=0,nBootInference=2000,useInfBootForTuning=T,alpha=0.05, TuningScheme="Optimal",OmegaOutlierThresh=10^(10),tauQuantReg=NULL){
  
  #Helpful variable definitions for bootstrap 
  m <- nrow(ProxyDat)
  p <- ncol(GoodDat)-1
 
  
  #Check formatting and print warning messages
  if(sum(PiLabel==0)>0){print("Warning: the probabilities that each sample 
                              is labelled should be strictly positive for statistical guarantees")}
  
  IdxCalib <- complete.cases(GoodDat)
  if(!is.null(clusterID)){
    ClustersCalib <- clusterID[IdxCalib]
    ClustersNonCalib <- clusterID[!IdxCalib]
    if(length(dplyr::intersect(ClustersCalib,ClustersNonCalib))>0){
      print("Warning: some clusters have a combination of calibration and noncalibration samples. Method assumes each cluster is either only in the either entirely in the calibration sample or not in the calibration sample at all.")
    }
    
    if(is.null(PiLabl)){
      print("Warning: PiLabel set to NULL, this can cause issues.") 
      print("If the labelled data is from the same distribution as the unlabelled a constant vector for PiLabel with values n/m where m is the totatl sample size and n is number of calibration samples (or the expected number of calibration samples when calibration sample size is randomly selected)")
    }
    
    if(!is.null(PiLabel)){
      LabelClusterSummaries <- data.frame(PiLabel,clusterID) %>% group_by(clusterID) %>% summarise(piLabMax=max(PiLabel),piLabMin=min(PiLabel),PiLabelCluster=mean(PiLabel))
        if(max(LabelClusterSummaries$piLabMax-LabelClusterSummaries$piLabMin)>0){
          print("Warning: label probabilities must be the same within each cluster")
        }
      }
    if((sum(BootstrapScheme %in% c("FullBootStudentized","CLTBased"))>0) ){
      print("Warning: Function doesn't currently support Studentized bootstrap, CLT based confidence intervals or CLT based tuning for clustered samples.")
    }
    if(RegType %in% c("Quantile","Quantile Regression")){
      print("Warning: Sandwich estimators for quantile regression with clustered samples not available. Some of the standard errors could be wrong")
    }
  }
  
  if(RegType %in% c("Quantile","Quantile Regression")){
    if((sum(BootstrapScheme %in% c("FullBootStudentized","CLTBased"))>0) ){
      print("Warning: Function doesn't support Studentized bootstrap, CLT based confidence intervals or CLT based tuning for quantile regression.")
    }
  } 
  
  
  
  #Extract Good and Proxy Data on the Complete cases 
  GoodDatCalib <- GoodDat[IdxCalib,]
  ProxyDatCalib <- ProxyDat[IdxCalib,]
  PiCalib <- PiLabel[IdxCalib]
  
  #Extract data on the noncomplete cases
  ProxyDatNonCalib <- ProxyDat[!IdxCalib,]
  PiNonCalib <- PiLabel[!IdxCalib]
  
  #Calculating Point estimates and standard errors
  
  
  #betaHat
  betaHatFit <- FormatAndFitReg(DatInp2 = GoodDatCalib,WeightsInp =  1/PiCalib,clusterIDInp2 = clusterID[IdxCalib],
                                RegTypeInp2 = RegType,OutcomeVarInp2 = OutcomeVarName,CalcSandwhich = T,quantileTauInp = tauQuantReg)
  betaHat <- betaHatFit$CoefsEst
  CoefNamesBetaHat <- names(betaHat)
  SandwichSEbetaHatCalib <- sqrt(diag(betaHatFit$VCOVSandwich))
  
  
  #GammaHat
  gammaHatCalibFit <- FormatAndFitReg(DatInp2 = ProxyDatCalib,WeightsInp  = 1/PiCalib,clusterIDInp2 = clusterID[IdxCalib],
                                      RegTypeInp2 = RegType,OutcomeVarInp2 = OutcomeVarName,CalcSandwhich = T,quantileTauInp = tauQuantReg)
  gammaHatCalib <- gammaHatCalibFit$CoefsEst
  SandwichCovgammaHatCalib <- gammaHatCalibFit$VCOVSandwich
  
  
  #gamma hat on rest of data 
  gammaHatNonCalibFit <- FormatAndFitReg(DatInp2 = ProxyDatNonCalib,WeightsInp =   1/(1-PiNonCalib),clusterIDInp2 = clusterID[!IdxCalib],
                                         RegTypeInp2 = RegType,OutcomeVarInp2 = OutcomeVarName,CalcSandwhich = T,quantileTauInp = tauQuantReg)
  gammaHatNonCalib <- gammaHatNonCalibFit$CoefsEst
  SandwichCovgammaHatNonCalib <- gammaHatNonCalibFit$VCOVSandwich
  
  
#########    Running Bootstrap ######################################
  
  #Collecting Bootstrap samples
  nBootTotal <- nBootTune+nBootInference
  if(sum(BootstrapScheme %in% c("FullBootStudentized"))>0){
    Boots <- BootstrapAndCalcEsts(ProxyDatInp = ProxyDat,GoodDatInp = GoodDat,PiLabelInp = PiLabel,clusterIDInp=clusterID,
                                      OutcomeVarInp = OutcomeVarName,RegTypeInp = RegType,B = nBootTotal,CalcNonCalib = TRUE,CalcPPSe=TRUE,quantTauInp2 = tauQuantReg)
  } else if(sum(BootstrapScheme %in% c("FullBootPercentile","FullBasicBoot"))>0){
    Boots <- BootstrapAndCalcEsts(ProxyDatInp = ProxyDat,GoodDatInp = GoodDat,PiLabelInp = PiLabel,clusterIDInp=clusterID,
                                      OutcomeVarInp = OutcomeVarName,RegTypeInp = RegType,B = nBootTotal,CalcNonCalib = TRUE,quantTauInp2 = tauQuantReg)
  } else if(sum(BootstrapScheme %in% c("QuickConvolution"))>0){
    Boots <- BootstrapAndCalcEsts(ProxyDatInp = ProxyDat,GoodDatInp = GoodDat,PiLabelInp = PiLabel,clusterIDInp=clusterID,
                                      OutcomeVarInp = OutcomeVarName,RegTypeInp = RegType,B = nBootTotal,CalcNonCalib = FALSE,quantTauInp2 = tauQuantReg)
  }
  if(useInfBootForTuning){
    idxBootTune <- 1:nBootTotal
    idxBootInf <- 1:nBootInference
  } else{
    idxBootTune <- 1:nBootTune
    idxBootInf <- (nBootTune+1):nBootTotal
  }
  
  ############ Constructing Confidence intervals for selected approaches ################################
  PPCIList <- list()
  OptOmegasList <- list()
  for(BootstrapSchemeCurr in BootstrapScheme){
    
    #Calculating relevant matrices for determining optimal tuning matrix
    if(BootstrapSchemeCurr=="CLTBased"){
      CalibBetaGammaNonBoot <- EstBetaGammaCalib_andVCOV.glm(GoodDatInp = GoodDatCalib,ProxyDatInp = ProxyDatCalib,PiInp = PiCalib,RegType = RegType,OutcomeVar = OutcomeVarName,Minp = m)
      A12hat <- CalibBetaGammaNonBoot$VCOVBetaGamma[1:(1+p),(2+p):(2*(1+p))]
      A22hat <- CalibBetaGammaNonBoot$VCOVBetaGamma[(2+p):(2*(1+p)),(2+p):(2*(1+p))]
      A33hat <- m*SandwichCovgammaHatNonCalib
    }
    
    #Calculating relevant matrices for determining optimal tuning matrix
    if((BootstrapSchemeCurr %in% c("FullBootStudentized","FullBootPercentile","FullBasicBoot")) & (sum(TuningScheme !="None")>0)){
      A12hat <- m*cov(Boots$betaHatCalibBoot[idxBootTune,],Boots$gammaHatCalibBoot[idxBootTune,])
      A22hat <- m*cov(Boots$gammaHatCalibBoot[idxBootTune,])
      A33hat <- m*cov(Boots$gammaHatNonCalibBoot[idxBootTune,])
    }
    
    if((BootstrapSchemeCurr %in% c("QuickConvolution")) & (sum(TuningScheme !="None")>0)){
      A12hat <- m*cov(Boots$betaHatCalibBoot[idxBootTune,],Boots$gammaHatCalibBoot[idxBootTune,])
      A22hat <- m*cov(Boots$gammaHatCalibBoot[idxBootTune,])
      A33hat <- m*SandwichCovgammaHatNonCalib
    }
    
    for(TuneSchemeCurr in TuningScheme){
     
      ##################### Calculating desired tuning matrices ####################
      if(TuneSchemeCurr=="None"){OmegaTuned <- diag(p+1)}
      if(TuneSchemeCurr=="Zero"){OmegaTuned <- matrix(0,nrow= (p+1),ncol= (p+1))}
      if(TuneSchemeCurr == "Diagonal"){OmegaTuned <- diag(diag(A12hat)/diag(A22hat+A33hat))}
      if(TuneSchemeCurr == "Optimal"){
        if(rcond(A22hat+A33hat) > (.Machine$double.eps)){
            OmegaTuned <- A12hat %*% solve(A22hat+A33hat)
        } else{
            OmegaTuned <- matrix(Inf,nrow = ncol(Boots$betaHatCalibBoot),ncol = ncol(Boots$gammaHatCalibBoot))
        }
      }
      
      #Making sure that OmegaOptEstimate didn't blow up and setting to zero if it did.
      if(mean(abs(OmegaTuned))>OmegaOutlierThresh){
        print("Warning: Optimal omega estimate below likely to large. Setting to zero.")
        print(paste0("Current tuning scheme:",TuneSchemeCurr,"Opt Omega Est:"))
        print(OmegaTuned)
        if(TuningScheme!="CLTBasedTune"){print("Consider increasing nBootTune")}
        dimOmegaOptEst <- dim(OmegaTuned)
        OmegaTuned <- matrix(0,nrow = dimOmegaOptEst[1],ncol = dimOmegaOptEst[2])
     }
      
      
      #Prediction-Powered point estimate
      betaPPest <- as.vector(betaHat+ OmegaTuned %*% (gammaHatNonCalib-gammaHatCalib))
      
      
      #Calculating asymptotic variance using CLT (if necessary)
      if(sum(BootstrapScheme %in% c("FullBootStudentized","CLTBased"))>0){
        CalibBetaGammaNonBoot <- EstBetaGammaCalib_andVCOV.glm(GoodDatInp = GoodDatCalib,ProxyDatInp = ProxyDatCalib,PiInp = PiCalib,RegType = RegType,OutcomeVar = OutcomeVarName,Minp = m)
        A11hat <- CalibBetaGammaNonBoot$VCOVBetaGamma[1:(1+p),1:(1+p)]
        A12hat <- CalibBetaGammaNonBoot$VCOVBetaGamma[1:(1+p),(2+p):(2*(1+p))]
        A22hat <- CalibBetaGammaNonBoot$VCOVBetaGamma[(2+p):(2*(1+p)),(2+p):(2*(1+p))]
        A33hat <- m*SandwichCovgammaHatNonCalib
        AsympCovPP <- OmegaTuned %*% (A22hat+A33hat) %*% t(OmegaTuned) - A12hat %*%  t(OmegaTuned)- t(A12hat %*%  t(OmegaTuned))+A11hat
        SEPPAsympEst <- sqrt(diag(AsympCovPP/m))
      }
      
      
      ##################### Calculating Confidence Interval for fixed tuning matrix choice ####################
      #CI that doesn't use bootstrap at all
      if(BootstrapSchemeCurr=="CLTBased"){
        CI_lbs <- matrix(betaPPest-qnorm(1-alpha/2)*SEPPAsympEst,ncol = 1)
        CI_ubs <- matrix(betaPPest+qnorm(1-alpha/2)*SEPPAsympEst,ncol = 1)
        CICurrScheme <- cbind(matrix(betaPPest,ncol=1),CI_lbs,CI_ubs)
      }
      
      #CI based on convolution bootstrap (Algorithm 3)
      if(BootstrapSchemeCurr=="QuickConvolution"){
        CalibDatBiasEstBoot <- Boots$betaHatCalibBoot[idxBootInf,] - t(OmegaTuned %*% t(Boots$gammaHatCalibBoot[idxBootInf,]))

        IIDGaussianNoise <- matrix(rnorm(nrow(CalibDatBiasEstBoot)*ncol(CalibDatBiasEstBoot)),nrow = ncol(CalibDatBiasEstBoot),ncol = nrow(CalibDatBiasEstBoot))
        AympNoiseNonCalibTerm <- t(OmegaTuned %*% t(chol(SandwichCovgammaHatNonCalib)) %*% IIDGaussianNoise)
        OmegaGammaHatNonCalibTerm <- matrix(rep(OmegaTuned %*% gammaHatNonCalib,nBootInference),nrow = nBootInference,byrow = T)
        QuickPPEstsBoot <- OmegaGammaHatNonCalibTerm + AympNoiseNonCalibTerm+CalibDatBiasEstBoot
        CILBSUBS <- t(sapply(data.frame(QuickPPEstsBoot), FUN = function(x) quantile(x,probs = c(alpha/2,1-alpha/2))))
        CICurrScheme <- cbind(matrix(betaPPest,ncol=1),CILBSUBS)
      }
      
      
      #CI based on percentile bootstrap (Algorithm 2)
      if(BootstrapSchemeCurr=="FullBootPercentile"){
        PPEstBoot <-  Boots$betaHatCalibBoot[idxBootInf,] + t(OmegaTuned %*% t(Boots$gammaHatNonCalibBoot[idxBootInf,]-Boots$gammaHatCalibBoot[idxBootInf,]))
        CILBSUBS <- t(sapply(data.frame(PPEstBoot), FUN = function(x) quantile(x,probs = c(alpha/2,1-alpha/2))))
        CICurrScheme <- cbind(matrix(betaPPest,ncol=1),CILBSUBS)
      }
      
      #CI based on basic bootstrap
      if(BootstrapSchemeCurr=="FullBasicBoot"){
        PPEstBoot <-  Boots$betaHatCalibBoot[idxBootInf,] + t(OmegaTuned %*% t(Boots$gammaHatNonCalibBoot[idxBootInf,]-Boots$gammaHatCalibBoot[idxBootInf,]))
        PivotsUnstd <- PPEstBoot-matrix(rep(betaPPest,nBootInference),nrow = nBootInference,byrow = T)
        PivotsUnstdQuants <- t(sapply(data.frame(PivotsUnstd), FUN = function(x) quantile(x,probs = c(alpha/2,1-alpha/2))))
        CICurrScheme <- cbind(matrix(betaPPest,ncol=1),betaPPest-PivotsUnstdQuants[,2],betaPPest-PivotsUnstdQuants[,1])
      }
      
      #CI based on studentized bootstrap
      if(BootstrapSchemeCurr=="FullBootStudentized"){
        PPEstBoot <- Boots$betaHatCalibBoot + t(OmegaTuned %*% t(Boots$gammaHatNonCalibBoot-Boots$gammaHatCalibBoot))
        PPSEBoot <- matrix(data = NA,nrow = nBootInference,ncol = (p+1))
        for(b in idxBootInf){
          AsympCovPPCurrBoot <- OmegaTuned %*% (Boots$ACovList[[b]]$A22hat+Boots$ACovList[[b]]$A33hat) %*% t(OmegaTuned) - Boots$ACovList[[b]]$A12hat %*%  t(OmegaTuned)- t(Boots$ACovList[[b]]$A12hat %*%  t(OmegaTuned))+ Boots$ACovList[[b]]$A11hat
          PPSEBoot[b,] <- sqrt(diag(AsympCovPPCurrBoot/m))
        }
        StandardizedPivots <- (PPEstBoot-matrix(rep(betaPPest,nBootInference),nrow = nBootInference,byrow = T))/PPSEBoot
        PivotQuantiles <- t(sapply(data.frame(StandardizedPivots), FUN = function(x) quantile(x,probs = c(alpha/2,1-alpha/2))))
        CICurrScheme <- cbind(matrix(betaPPest,ncol=1),betaPPest-SEPPAsympEst*PivotQuantiles[,2],betaPPest-SEPPAsympEst*PivotQuantiles[,1])
      }
      

      
      #################### Formatting Output ##########################
      rownames(CICurrScheme) <- CoefNamesBetaHat
      colnames(CICurrScheme) <- c("Estimate","CI_lb","CI_ub") 
      PPCIList[[paste0(BootstrapSchemeCurr,"_",TuneSchemeCurr)]] <- CICurrScheme
      OptOmegasList[[paste0(BootstrapSchemeCurr,"_",TuneSchemeCurr)]] <- OmegaTuned
      CICurrScheme <- NULL
      OmegaTuned <- NULL

    }
  }
  
  #Alternative beta estimates and CIs
  
  #naive approach
  
  gammaHatAllFit <- FormatAndFitReg(DatInp2 = ProxyDat,RegTypeInp2 = RegType,OutcomeVarInp2 = OutcomeVarName,CalcSandwhich = T,quantileTauInp = tauQuantReg)
  gammaHatAll <- gammaHatAllFit$CoefsEst
  SandwichSEgammaHatAll <- sqrt(diag(gammaHatAllFit$VCOVSandwich))
  
  GammaHatEstCI <- cbind(gammaHatAll,
                         gammaHatAll-qnorm(1-alpha/2)*as.vector(SandwichSEgammaHatAll),
                         gammaHatAll+qnorm(1-alpha/2)*as.vector(SandwichSEgammaHatAll))
  rownames(GammaHatEstCI) <- names(gammaHatAll)
  colnames(GammaHatEstCI) <- c("Estimate","CI_lb","CI_ub")
  PPCIList[["naive"]] <- GammaHatEstCI

  
  #classical approach
  betaHatEstCI <- cbind(betaHat,
                         betaHat-qnorm(1-alpha/2)*as.vector(SandwichSEbetaHatCalib),
                         betaHat+qnorm(1-alpha/2)*as.vector(SandwichSEbetaHatCalib))
  rownames(betaHatEstCI) <- names(betaHat)
  colnames(betaHatEstCI) <- c("Estimate","CI_lb","CI_ub")
  PPCIList[["classical"]] <- betaHatEstCI

  return(list(OptOmegasList=OptOmegasList,CIList=PPCIList))
}