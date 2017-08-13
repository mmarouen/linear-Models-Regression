extractLambdas <- function(X,DoF=NULL){
  s=svd(X)
  dj=s$d
  p=ncol(X)
  if(! is.null(DoF)){range=DoF}
  else(range=seq(0,p,0.25))
  N=length(range)
  Lambdas=vector(length = N)
  Lambda1=0#initialize lambda=0
  for (i in 1:N){
    k=range[N-i+1]
    Lambda=Lambda1
    d_lambda=sum(dj^2/(dj^2+Lambda))-k
    d_lambda_prime=-sum(dj^2/((dj^2+Lambda)^2))
    Lambda1=Lambda-d_lambda/d_lambda_prime
    while(abs(d_lambda)>0.0001){
      Lambda=Lambda1
      d_lambda=sum(dj^2/(dj^2+Lambda))-k
      d_lambda_prime=-sum(dj^2/((dj^2+Lambda)^2))
      Lambda1=Lambda-d_lambda/d_lambda_prime
    }
    Lambdas[i]=Lambda1
  }
  return(Lambdas)
}

#Select best subset with regards a loss criterion by cross validation
SubsetSelect <- function(Input,Response,size=NULL,CV=FALSE,cvScale=TRUE){
  
  #init parameters:
  errorVector=c()#CV errors
  betaMat=c()#in case no DoF or shrinkage will return a matrix of all coefficients else empty
  beta=c()#in case no DoF or shrinkage will return beta coefficient else will return beta which give lowest error
  error=c()#RSS related to beta
  yhat=c()#training prediction
  varList=c()#in case no DoF or shrinkage empty, else will return variables which yield lowest RSS
  yhat0=c()
  modelF=FALSE#if TRUE, returns the best model
  stdDev=c()#standard deviations for CV
  DoF=c()

  if(!is.null(size)){q=floor(size)}
  if(is.null(size)){
    q=ncol(Input)
    modelF=TRUE
  }
  if(!CV){
    Input=as.matrix(Input)
    N=nrow(Input)
    p=ncol(Input)
    DoF=0:p
    out=subsetfit(Input,Response,q=q,modelFit = modelF)
    beta=out$beta
    RSS=out$RSS
    betaMat=out$betaMat
    errorVector=out$errorVector
    yhat=predictLS(X=Input,beta=beta)
  }
  if(CV){
    source("D:/RProject/toolkits/ModelAssessment&Selection.R")
    Input=as.matrix(Input)
    N=nrow(Input)
    p=ncol(Input)
    DoF=0:p
    for (k in 1:(p+1)){
      out=crossValidate(Input,Response,type="Subset",scaleInputs=cvScale,scaleType="Standardize",complexity=DoF[k])
      errVector=out$errorVector
      stdDev[k]=out$stdDev
      errorVector[k]=mean(errVector)
    }
    j0=which.min(errorVector)-1
    if(!cvScale){
      Input=ScaleData(Input,Response,type = "Standardize")$train
    }
    out00=subsetfit(Input,Response,q = j0)
    beta=out00$beta
    yhat=predictLS(X=Input,beta=beta)
    error=mean((yhat-Response)^2)
    varList=which(beta[-1]!=0)
  }

  return(list(yhat=yhat,beta=beta,RSS=error,betaMat=betaMat,varList=varList,
              errorVector=errorVector,stdDev=stdDev,DoF=DoF))
}

#Select best subset with a given size constraint
subsetfit<-function(Input,Response,q,modelFit=FALSE){
  
  #init parameters
  Input=as.matrix(Input)
  p=ncol(Input)
  N=nrow(Input)
  betaMat=matrix(0,nrow=p+1);betaMat[1,]=mean(Response)
  beta=rep(0,length=p+1);beta[1]=mean(Response)
  RSS=0
  errorVector=rep(0,length=q+1)
  yhat=0
  featList=0
  
  if(q==0){
    yhat=mean(Response)
    RSS=mean((yhat-Response)^2)
  }
  if(q>0){
    count=q
    if(modelFit){
      count=1
      errorVector=rep(0,length=q+1)
      errorVector[1]=mean((Response-mean(Response))^2)
    }
    for(k in count:q){
      set=combn(p,k)
      Nk=ncol(set)
      error=rep(1,length=Nk)
      beta0=matrix(0,ncol=Nk,nrow=k)
      yhat=mean(Response)
      for(j in 1:Nk){
        featList0=set[,j]
        X=Input[,featList0]
        beta0[,j]=solve(t(X)%*%X)%*%t(X)%*%(Response-mean(Response))
        if(k==1){yhat0=beta0[,j]*X+mean(Response)}
        if(k>1){yhat0=X%*%beta0[,j]+mean(Response)}
        error[j]=mean((yhat0-Response)^2)
      }
      j0=which.min(error)
      featList=set[,j0]
      RSS=error[j0]
      beta[featList+1]=beta0[,j0]
      if(modelFit){
        errorVector[k+1]=RSS
        betaMat=cbind(betaMat,beta)
      }
    }
    if(modelFit){
      j0=which.min(errorVector)
      RSS=min(errorVector)
      beta=betaMat[,j0]
      featList=which(beta[-1]!=0)
    }
    yhat=cbind(1,Input)%*%beta
  }
  return(list(beta=beta,RSS=RSS,yhat=yhat,varList=featList,errorVector=errorVector,betaMat=betaMat))
}
