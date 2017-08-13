#select best PLS by CV
PLS <- function(Input,resp,df=NULL,CV=FALSE,cvScale=TRUE){
  
  #init parameters:
  errVector=c()#CV errors
  betaMat=c()#in case no DoF or shrinkage will return a matrix of all coefficients else empty
  beta=c()#in case no DoF or shrinkage will return beta coefficient else will return beta which give lowest error
  error=c()#RSS related to beta
  yhat=c()#training prediction
  varList=c()#in case no DoF or shrinkage empty, else will return variables which yield lowest RSS
  yhat0=c()
  stdDev=c()#standard deviations for CV
  DoF=c()
  yhatMat=c()
  
  if(!is.null(df)){q=min(floor(df),ncol(Input))}
  if(is.null(df)){
    q=ncol(Input)
    modelF=TRUE
  }
  if(!CV){
    out=PLSfit(Input,resp,q)
    beta=out$Beta
    RSS=out$RSS
    betaMat=out$betaMat
    errorVector=out$errorVector
    yhatMat=out$yhatMat
    yhat=predictLS(X=Input,beta=beta)
    DoF=0:q
  }
  if(CV){
    source("D:/RProject/toolkits/ModelAssessment&Selection.R")
    DoF=0:q
    out=crossValidate(Input,resp,type="PLS",scaleInputs=cvScale,
                      scaleType="Standardize",complexity = DoF)
    errorVector=out$errorVector
    stdDev=out$stdDevVector
    error=min(errorVector)
    bestSize=which.min(error)-1
    betaMat=out$betaMat
    if(cvScale){Input=ScaleData(Input,resp,type = "Standardize")$train}
    out00=PLSfit(Input,resp,bestSize)
    beta=out00$Beta
    yhat=predictLS(X=Input,beta=beta)
  }
  return(list(yhat=yhat,beta=beta,yhatMat=yhatMat,RSS=error,betaMat=betaMat,errorVector=errorVector,stdDev=stdDev,DoF=DoF))
}

#PLS
PLSfit<- function(Input,resp,q){
  
  resp0=resp-mean(resp)
  Input=as.matrix(Input)
  p=ncol(Input)
  N=nrow(Input)
  betaMat=matrix(0,ncol=q+1,nrow=p+1)
  betaMat[1,]=mean(resp)
  errorVector=rep(0,length=q+1)
  yhatMat=matrix(0,nrow=N,ncol=q+1)
  yhatMat[,1]=mean(resp)
  Beta=c()
  RSS=c()
  km=matrix()
  Km=diag(p)
  
  if(q==0){
    Beta=betaMat[,1]
    RSS=mean((resp-mean(resp))^2)
  }
  if(q>0){
    Xj=Input
    RSS[1]=mean((resp-mean(resp))^2)
    for(j in 1:q){
      u=t(Xj)%*%resp0
      z=Xj%*%u
      Z2=sum(z^2)
      theta=sum(z*resp)/Z2
      yhatMat[,j+1]=yhatMat[,j]+theta*z
      betaMat[2:(p+1),j+1]=betaMat[2:(p+1),j]+theta*Km%*%Km%*%t(Xj)%*%resp0
      km=diag(p)-u%*%t(u)%*%t(Xj)%*%Xj/Z2
      Xj=Xj%*%km
      Km=Km%*%km
      RSS[j+1]=mean((resp-yhatMat[,j+1])^2)
    }
    Beta=betaMat[,which.min(RSS)]
  }
  return(list(yhatMat=yhatMat,betaMat=betaMat,RSS=RSS,Beta=Beta))
}

#predictions for any linear model
predictLS<-function(LSModel=NULL,X,beta=NULL){
  if(!is.null(LSModel)){beta=LSModel$beta}
  yhat=as.matrix(cbind(1,X))%*%beta
  return(yhat)
}

