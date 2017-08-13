#Ridge selection
Ridge <-function(Input,resp,Lambdas=NULL,df=NULL,CV=FALSE,cvScale=TRUE){
  
  #init parameters
  errorVector=c()
  RSS=c()
  yhat=c()
  beta=as.matrix(rep(0,length=ncol(Input)+1))
  betaMat=matrix(0,nrow = ncol(Input)+1)
  rownames(betaMat)=c("intercept",colnames(Input))
  stdDev=c()
  DoF=c()
  
  if(!is.null(Lambdas) | !is.null(df)){
    if(!is.null(Lambdas)){ll=Lambdas}
    if(!is.null(df) & is.null(Lambdas)){ll=extractLambdas(Input,DoF = df)}
    out=ridgefit(Input,resp,lambda = ll)
    RSS=out$RSS
    beta=out$beta
    yhat=out$yhat
  }
  
  if(is.null(df) & is.null(Lambdas)){
    lambda_vector=extractLambdas(Input)
    L=length(lambda_vector)
    N=nrow(Input)
    p=ncol(Input)
    Input=as.matrix(Input)
    for (i in 1:L){
      Lambda=lambda_vector[i]
      if(CV){
        source("D:/RProject/toolkits/ModelAssessment&Selection.R")
        out=crossValidate(Input,resp,type = "Ridge",complexity = Lambda,scaleInputs = cvScale,scaleType ="Standardize")
        errorVector[i]=mean(out$errorVector)
        stdDev[i]=out$stdDev
      }
      out0=ridgefit(Input,resp,Lambda)
      if(!CV){errorVector[i]=out0$RSS}
      betaMat=cbind(betaMat,out0$beta)
      Hat=Input%*% solve(t(Input)%*%Input + Lambda * diag(p)) %*% t(Input)
      DoF[i]=sum(diag(Hat))
    }
    betaMat=betaMat[,-1]
    beta=betaMat[,which.min(errorVector)]
    RSS=min(errorVector)
    yhat=predictLS(X = Input,beta = beta)
  }
  
  return(list(yhat=yhat,beta=beta,RSS=RSS,betaMat=betaMat,errorVector=errorVector,
              stdDev=stdDev,DoF=DoF))
}

#Builds a ridge model given lambda or df
ridgefit <-function(Input,resp,lambda){
  X=as.matrix(Input)
  resp0=resp-mean(resp)
  mat=solve(t(X)%*%X+lambda*diag(ncol(X)))
  beta=c(mean(resp),mat%*%t(X)%*%resp0)
  yhat=as.matrix(cbind(1,X))%*%beta
  RSS=mean((resp-yhat)^2)
  return(list(yhat=yhat,beta=beta,RSS=RSS))
}

#converts degree of freedom into L2 penalty coefficient relatively to a dataset
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
