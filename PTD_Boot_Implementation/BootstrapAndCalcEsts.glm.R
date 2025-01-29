#Helper functions


#Resample B times and calculate hatGammaCalib and hatBetaCalib on each resampling 
#(if CalcNonCalib=TRUE also calculate hatGammaNonCalib)
BootstrapAndCalcEsts.glm <- function(ProxyDatInp,GoodDatInp,PiLabelInp,clusterIDInp,OutcomeVarInp,RegTypeInp,B,CalcNonCalib=FALSE,CalcPPSe=FALSE,quantTauInp2=NULL){
  m <- nrow(ProxyDatInp)
  p <- ncol(GoodDatInp)-1
  
  
  #Formatting Clusters for cluster sampling
  if(!is.null(clusterIDInp)){
    uniqueClusterIDs <- unique(clusterIDInp)
    idxInEachCluster <- list()
    for(i in 1:length(uniqueClusterIDs)){
      idxInEachCluster[[i]] <- which(clusterIDInp==uniqueClusterIDs[i])
    }
  }
  
  
  gammaHatCalibBoot <- matrix(data = NA,nrow = B,ncol = (p+1))
  betaHatCalibBoot <-  matrix(data = NA,nrow = B,ncol = (p+1))
  gammaHatNonCalibBoot <- matrix(data = NA,nrow = B,ncol = (p+1))
  ACovList <- list()
  
  
  
  
  if (CalcNonCalib){
    #Bootstrap full data (Slower)
    for(i in 1:B){
      
      if(!is.null(clusterIDInp)){ #Clustered bootstrap
        IdxClusterBoot <- sample(1:length(uniqueClusterIDs),size = length(uniqueClusterIDs),replace = T)
        idxBoot <- vector()
        for(j in 1:length(IdxClusterBoot)){
          idxBoot <- c(idxBoot,idxInEachCluster[[IdxClusterBoot[j]]])
        }
      }else{ #Regular bootstrap
        idxBoot <- sample(x = 1:m,size = m,replace = T)
      }
      
      #Defining Bootstrapped data
      ProxyDatBoot <- ProxyDatInp[idxBoot,]
      GoodDatBoot <- GoodDatInp[idxBoot,]
      PiLabelBoot <- PiLabelInp[idxBoot]
      IdxCalibBoot <- complete.cases(GoodDatBoot)
      
      betaHatCalibBoot[i,] <- FormatAndFitReg(DatInp2 = GoodDatBoot[IdxCalibBoot,],WeightsInp =  1/PiLabelBoot[IdxCalibBoot],RegTypeInp2 = RegTypeInp,
                                              OutcomeVarInp2 = OutcomeVarInp,CalcSandwhich = F,quantileTauInp = quantTauInp2)$CoefsEst
      gammaHatCalibBoot[i,] <- FormatAndFitReg(DatInp2 = ProxyDatBoot[IdxCalibBoot,],WeightsInp =  1/PiLabelBoot[IdxCalibBoot],RegTypeInp2 = RegTypeInp,
                                               OutcomeVarInp2 = OutcomeVarInp,CalcSandwhich = F,quantileTauInp = quantTauInp2)$CoefsEst
      NonCalibBootFitCurr <- FormatAndFitReg(DatInp2 = ProxyDatBoot[!IdxCalibBoot,],WeightsInp =  1/(1-PiLabelBoot[!IdxCalibBoot]),RegTypeInp2 = RegTypeInp,
                      OutcomeVarInp2 = OutcomeVarInp,CalcSandwhich = CalcPPSe,quantileTauInp = quantTauInp2)
      gammaHatNonCalibBoot[i,] <- NonCalibBootFitCurr$CoefsEst

      if(CalcPPSe){ #Only use for studentized bootstrap
        ACovList[[i]] <- list()
        CalibBetaGammaBoot <- EstBetaGammaCalib_andVCOV.glm(GoodDatInp = GoodDatBoot[IdxCalibBoot,],ProxyDatInp = ProxyDatBoot[IdxCalibBoot,],PiInp = PiLabelBoot[IdxCalibBoot],RegType = RegTypeInp,OutcomeVar = OutcomeVarInp,Minp = m)
        ACovList[[i]]$A11hat <- CalibBetaGammaBoot$VCOVBetaGamma[1:(1+p),1:(1+p)]
        ACovList[[i]]$A12hat <- CalibBetaGammaBoot$VCOVBetaGamma[1:(1+p),(2+p):(2*(1+p))]
        ACovList[[i]]$A22hat <- CalibBetaGammaBoot$VCOVBetaGamma[(2+p):(2*(1+p)),(2+p):(2*(1+p))]
        ACovList[[i]]$A33hat <- m*NonCalibBootFitCurr$VCOVSandwich
      } 
    }
    
  } else if (CalcNonCalib==FALSE){
    #Bootstrap just the calibration sample (quicker)
    IdxCalibInp <- complete.cases(GoodDatInp)
    GoodDatCalibInp <- GoodDatInp[IdxCalibInp,]
    ProxyDatCalibInp <- ProxyDatInp[IdxCalibInp,]
    PiCalibInp <- PiLabelInp[IdxCalibInp]
    
    #Note we use Quicker way of bootstrapping calibration data that is equivalent 
    #to resampling all data (or clusters) with replacement and only selecting labelled samples among them
    if(!is.null(clusterIDInp)){
      ClusterCalibInp <- clusterIDInp[IdxCalibInp]
      uniqueClusterIDsCalib <- unique(ClusterCalibInp)
      idxInEachCalibCluster <- list()
      for(i in 1:length(uniqueClusterIDsCalib)){
        idxInEachCalibCluster[[i]] <- which(ClusterCalibInp==uniqueClusterIDsCalib[i])
      }
      nClustersCalibBoot <-  rbinom(n = B,size = length(uniqueClusterIDs),prob = length(uniqueClusterIDsCalib)/length(uniqueClusterIDs))
      #print(table(nClustersCalibBoot))
    } else{
      nCalibBoot <- rbinom(n = B,size = m,prob = nrow(GoodDatCalibInp)/m)
      
    }
    
    
    for(i in 1:B){
      if(!is.null(clusterIDInp)){ #clustered bootstrap
        IdxClusterCalibBoot <- sample(1:length(uniqueClusterIDsCalib),size = nClustersCalibBoot[i],replace = T)
        idxBootForCalib <- vector()
        for(j in 1:length(IdxClusterCalibBoot)){
          idxBootForCalib <- c(idxBootForCalib,idxInEachCalibCluster[[IdxClusterCalibBoot[j]]])
        }
      }else{ #regular bootstrap
        idxBootForCalib <- sample(x = 1:nrow(GoodDatCalibInp),size = nCalibBoot[i],replace = T)
      }
      
      GoodDatCalibBoot <- GoodDatCalibInp[idxBootForCalib,]
      ProxyDatCalibBoot <- ProxyDatCalibInp[idxBootForCalib,]
      PiCalibBoot <- PiCalibInp[idxBootForCalib]
      betaHatCalibBoot[i,] <- FormatAndFitReg(DatInp2 = GoodDatCalibBoot,WeightsInp =  1/PiCalibBoot,RegTypeInp2 = RegTypeInp,
                                              OutcomeVarInp2 = OutcomeVarInp,CalcSandwhich = F,quantileTauInp = quantTauInp2)$CoefsEst
      gammaHatCalibBoot[i,] <- FormatAndFitReg(DatInp2 = ProxyDatCalibBoot,WeightsInp = 1/PiCalibBoot,RegTypeInp2 = RegTypeInp,
                                               OutcomeVarInp2 = OutcomeVarInp,CalcSandwhich = F,quantileTauInp = quantTauInp2)$CoefsEst
    }
  }
  
  return(list(betaHatCalibBoot=betaHatCalibBoot,gammaHatCalibBoot=gammaHatCalibBoot,gammaHatNonCalibBoot=gammaHatNonCalibBoot,ACovList=ACovList))
}
