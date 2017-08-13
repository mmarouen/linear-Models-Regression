#selects PCA reduction by CV
PCR <- function(Input,Response,df=NULL,CV=FALSE,cvScale=TRUE){
  
  #init parameters:
  errVector=c()#CV errors
  betaMat=c()#in case no DoF or shrinkage will return a matrix of all coefficients else empty
  beta=c()#in case no DoF or shrinkage will return beta coefficient else will return beta which give lowest error
  error=c()#RSS related to beta
  yhat=c()#training prediction
  varList=c()#in case no DoF or shrinkage empty, else will return variables which yield lowest RSS
  yhat0=c()
  modelF=FALSE#if TRUE, will the best model
  stdDev=c()#standard deviations for CV
  DoF=c()
  
  if(!is.null(df)){q=floor(df)}
  if(is.null(df)){
    q=ncol(Input)
    modelF=TRUE
  }
  if(!CV){
    out=PCRfit(Input,Response,q,modelFit = modelF)
    beta=out$beta
    RSS=out$RSS
    betaMat=out$betaMat
    errorVector=out$errorVector
    yhat=predictLS(X=Input,beta=beta)
    DoF=0:q
  }
  if(CV){
    source("D:/RProject/toolkits/ModelAssessment&Selection.R")
    DoF=0:ncol(Input)
    out=crossValidate(Input,Response,type="PCR",scaleInputs=cvScale,scaleType="Standardize",complexity = DoF)
    errorVector=out$errorVector
    stdDev=out$stdDevVector
    error=min(errorVector)
    bestSize=which.min(error)-1
    if(cvScale){Input=ScaleData(Input,Response,type = "Standardize")$train}
    out00=PCRfit(Input,Response,bestSize)
    beta=out00$beta
    yhat=predictLS(X=Input,beta=beta)
  }
  return(list(yhat=yhat,beta=beta,RSS=error,betaMat=betaMat,errorVector=errorVector,
              stdDev=stdDev,DoF=DoF))
}

#prediction function for PCR. Prediction can be done 02 ways:
#provide the best selected model OR training data + subset size
PCRfit <- function(Input,resp,q,modelFit=FALSE){
  
  resp0=resp-mean(resp)
  Input=as.matrix(Input)
  p=ncol(Input)
  N=nrow(Input)
  betaMat=matrix(0,ncol=q+1,nrow=p+1)
  betaMat[1,]=mean(resp)
  errorVector=rep(0,length=q+1)
  RSS=c()
  beta=c()
  yhat=c()
  
  pca=PCA(Input)
  Z=pca$Z;V=pca$V;D=pca$D
  if(q==0){
    RSS=mean((resp-mean(resp))^2)
    beta=betaMat[,1]
  }
  if(q>0){
    thetaHat=rep(0,length=q)
    if(modelFit){errorVector[1]=mean((resp-mean(resp))^2)}
    for (j in 1:q){
      thetaHat[j]=sum(t(Z[,j])%*%resp0)/D[j,j]
      for (i in j:q){
        betaMat[2:(p+1),i+1]=betaMat[2:(p+1),i+1]+thetaHat[j]*V[,j]
      }
      if(modelFit){
        yhat=predictLS(X=Input,beta = betaMat[,i+1])
        errorVector[j+1]=mean((yhat-resp)^2)
      }
    }
    if(modelFit){j0=which.min(errorVector)}
    if(!modelFit){j0=q+1}
    beta=betaMat[,j0]
    RSS=errorVector[j0]
    yhat=predictLS(X = Input,beta = beta)
  }
  return(list(beta=beta,RSS=RSS,yhat=yhat,errorVector=errorVector,betaMat=betaMat))
}

#predictions for any linear model
predictLS<-function(LSModel=NULL,X,beta=NULL){
  if(!is.null(LSModel)){beta=LSModel$beta}
  yhat=as.matrix(cbind(1,X))%*%beta
  return(yhat)
}
