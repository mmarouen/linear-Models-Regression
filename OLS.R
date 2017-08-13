LR<-function(input,response){
  X=input
  Y=response
  yhat=c()
  X=cbind(1,X)
  X=as.matrix(X)
  beta=solve(t(X)%*%X)%*%t(X)%*%Y
  yhat=X%*%beta
  sigma2=(1/(ncol(X)))*sum((response-yhat)^2)
  var_beta=solve(t(X)%*%X)*sigma2
  return(list(yhat=yhat,beta=beta,sigma2=sigma2,varBeta=var_beta))
}

#predictions for any linear model
predictLS<-function(LSModel=NULL,X,beta=NULL){
  if(!is.null(LSModel)){beta=LSModel$beta}
  yhat=as.matrix(cbind(1,X))%*%beta
  return(yhat)
}
