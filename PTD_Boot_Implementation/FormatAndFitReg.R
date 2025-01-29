#Helper functions

library(quantreg)

#Fit GLM or quantile regression
FormatAndFitReg <- function(DatInp2,RegTypeInp2,WeightsInp=NULL,clusterIDInp2=NULL,OutcomeVarInp2="Y",CalcSandwhich=F,quantileTauInp=NULL){

  VCOVSandwich <- NULL
  RegFormulaLocal <- as.formula(paste0(OutcomeVarInp2,"~ ."))
  
  if(RegTypeInp2 %in% c("Quantile","Quantile Regression")){#Quantile regression
    QuantileRegFit <- rq(formula = RegFormulaLocal,data = DatInp2,tau = quantileTauInp,weights = WeightsInp)
    CoefsEst <- QuantileRegFit$coefficients
    if(CalcSandwhich){VCOVSandwich <- summary(QuantileRegFit,se="ker",covariance = TRUE)$cov}
    
  } else{ #GLMs
    if(RegTypeInp2 %in% c("linear","Linear","Gaussian","gaussian","Normal","normal","OLS","ols")){
      familyLocal <- gaussian(link = "identity")
    } else if(RegTypeInp2 %in% c("logistic","Logistic","logit")){
      familyLocal <- binomial(link = "logit")
    } else if(RegTypeInp2 %in% c("Poisson","poisson")){
      familyLocal <- poisson(link = "log")
    }

    RegFit <- glm(formula = RegFormulaLocal,family = familyLocal,data = DatInp2,weights = WeightsInp)
    CoefsEst <- RegFit$coefficients
    if(CalcSandwhich){
      if(is.null(clusterIDInp2)){
        VCOVSandwich <- sandwich(RegFit)
      } else{
        VCOVSandwich <- vcovCL(RegFit,cluster = clusterIDInp2,sandwich = T)
      }
    }
  }
 
  

  return(list(CoefsEst=CoefsEst,VCOVSandwich=VCOVSandwich))
}