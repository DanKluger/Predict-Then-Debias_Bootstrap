#This helper function applies to GLMs only. It does not apply to clustered or weighted sampling settings either.
#It takes a dataset as input and computes a sandwich estimator for the asymptotic covariance matrix of (hatBetaCalib,hatGammaCalib)  
#It is called when using CLT-based confidence intervals in PTDBootModularized. It is also called in each bootstrap interation when using the studentized bootstrap

EstBetaGammaCalib_andVCOV.glm <- function(GoodDatInp,ProxyDatInp,PiInp,RegType="linear",OutcomeVar="Y",Ninp){
  N <- Ninp
  p <- ncol(GoodDatInp)-1
  
  piVecVal <- PiInp
  X_augVal <- cbind(rep(1,nrow(GoodDatInp)),as.matrix(GoodDatInp)[,which(names(GoodDatInp) != OutcomeVar )] )
  tildeX_augVal <- cbind(rep(1,nrow(ProxyDatInp)),as.matrix(ProxyDatInp)[,which(names(ProxyDatInp) != OutcomeVar )] )

  if(RegType %in% c("linear","Linear","Gaussian","gaussian","Normal","normal","OLS","ols")){
    psidot <- function(s){return(s)}
    psidotdot <- function(s){return(rep(1,length(s)))}
    familyInp <- gaussian(link = "identity")
  } else if(RegType %in% c("logistic","Logistic","logit")){
    psidot <- function(s){return(1/(1+exp(-s)))}
    psidotdot <- function(s){return(exp(s)/((1+exp(s))^2))}
    familyInp <- binomial(link = "logit")
  } else if(RegType %in% c("Poisson","poisson")){
    psidot <- function(s){return(exp(s))}
    psidotdot <- function(s){return(exp(s))}
    familyInp <- poisson(link = "log")
  }
  
  
  #Fitting relevant GLMs
  #Formatting Formula for regression
  RegFormulaInp <- as.formula(paste0(OutcomeVar,"~ ."))
  
 
  GlmFitGoodDatBoot <- glm(formula = RegFormulaInp,family = familyInp,data = GoodDatInp,weights = 1/piVecVal)
  GlmFitProxyDatValSampBoot <- glm(formula = RegFormulaInp,family = familyInp,data = ProxyDatInp,weights = 1/piVecVal)
  
  #Coefficients 
  betaHatBoot <- GlmFitGoodDatBoot$coefficients
  gammaHatValBoot <- GlmFitProxyDatValSampBoot$coefficients

  #Estimating asymptotic covariance submatrices
  # (with the exception of C12, these can also be computed 
  #using the bread and meat functions in the sandwhich package)
  
  
  #Now calculating D matrices (expected derivatives of estimating equation)
  #Letting W=I/pi where I is an indactor of being in the validation set
  #D11=E[W psidotdot(\beta_*^\tran X) XX^\tran]
  #D22= E[psidotdot(\gamma_*^\tran \tilde{X}) \tilde{X} \tilde{X}^\tran]
  
  
  weightsX <- psidotdot(X_augVal %*% betaHatBoot)/piVecVal
  weightsXTilde <- psidotdot(tildeX_augVal %*% gammaHatValBoot)/piVecVal
  
  D11hat <- t(matrix(rep(weightsX,p+1), ncol = (p+1),byrow = F)*X_augVal) %*% X_augVal/N
  D22hat <- t(matrix(rep(weightsXTilde,p+1), ncol = (p+1),byrow = F)*tildeX_augVal) %*% tildeX_augVal/N
  
  #C=cov(W* (psidot(\beta_*^\tran X)-Y)*X,
  #   W*(psidot(\gamma_*^\tran \tilde{X})-\tilde{Y})*\tilde{X})
  
  WeightedScoreVecXmatVal <- matrix(rep((GoodDatInp[[OutcomeVar]]-psidot( X_augVal %*% betaHatBoot))/piVecVal
                                        ,p+1),ncol = (p+1),byrow = F)*X_augVal
  
  
  C11hat <- (t(WeightedScoreVecXmatVal) %*% WeightedScoreVecXmatVal)/N
  
  WeightedScoreVecTildeXCalib<-  matrix(rep((ProxyDatInp[[OutcomeVar]]-psidot( tildeX_augVal %*% gammaHatValBoot))/piVecVal
                                            ,p+1),ncol = (p+1),byrow = F)*tildeX_augVal
  
  C22hat <- (t(WeightedScoreVecTildeXCalib) %*% WeightedScoreVecTildeXCalib)/N
  

  C12hat <- (t(WeightedScoreVecXmatVal) %*% WeightedScoreVecTildeXCalib)/N
  
  #Calculating Asymptotic covariance matrix for (hatBetaCalibBoot, hatGammaCalibBoot)
  D11Inv <- solve(D11hat)
  D22Inv <- solve(D22hat)
  VCOV11 <- D11Inv %*% C11hat %*% D11Inv
  VCOV12 <- D11Inv %*% C12hat %*% D22Inv
  VCOV21 <- D22Inv %*% t(C12hat) %*% D11Inv
  VCOV22 <- D22Inv %*% C22hat %*% D22Inv
  VCOVOut <- rbind(cbind(VCOV11,VCOV12),cbind(VCOV21,VCOV22))
  
  
  
  
  return(list(betaHatBoot=betaHatBoot,GammaHatCalibBoot=gammaHatValBoot,VCOVBetaGamma=VCOVOut)) 
}