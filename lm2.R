lm2<-function(input,response,fit="OLS",df=NULL,CV=FALSE,cvScale=FALSE,lambda=NULL){
  #init
  yhat=c()
  beta=c()
  betaMat=c()
  yhatMat=c()
  RSS=c()
  errorVector=c()
  DoF=c()
  stdDev=c()
  sigma2=c() #linear regression
  varBeta=c() #linear regression
  
  if(fit=="OLS"){out=LR(input,response)}
  if(fit%in%c("LAR","Lasso","IFSR")){out=LAR(input,response,type=fit,CV=CV,df=df,
                                             cvScale=cvScale)}
  if(fit=="Subset"){out=SubsetSelect(input,response,size=df,CV=CV,cvScale=cvScale)}
  if(fit=="Ridge"){out=Ridge(input,response,Lambdas=lambda,df=df,CV=CV,cvScale=cvScale)}
  if(fit=="PCR"){out=PCR(input,response,df=df,CV=CV,cvScale=cvScale)}
  if(fit=="PLS"){out=PLS(input,response,df=df,CV=CV,cvScale=cvScale)}
  yhat=out$yhat
  beta=out$beta
  if(fit != "OLS"){
    yhatMat=out$yhatMat
    betaMat=out$betaMat
    errorVector=out$errorVector
    stdDev=out$stdDev
    DoF=out$DoF
    RSS=out$RSS
    l1=list(yhat=yhat,beta=beta,RSS=RSS,yhatMat=yhatMat,betaMat=betaMat,
                errorVector=errorVector,stdDev=stdDev,DoF=DoF)
    if(fit%in%c("Subset","LAR")){
      l1=c(l1,varList=out$varList)
    }
  }
  if(fit=="OLS"){
    sigma2=out$sigma2
    varBeta=out$varBeta
    l1=list(yhat=yhat,beta=beta,sigma2=sigma2,varBeta=varBeta)
  }
  return(l1)
}
