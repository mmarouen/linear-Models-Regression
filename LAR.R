#3 model selection algorithms: 
#Least Angle Regression, Incremental Forward Stagewise Regression, Lasso
#If no cross validation, must centralize & unit norm for covariates
#If cross validation: data scaling will be done in the CV loop
LAR<-function(Input,resp,type="LAR",df=NULL,CV=FALSE,cvScale=TRUE){
  
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
  
  if(!is.null(df)){size=min(floor(df),ncol(Input))}
  if(is.null(df)){
    size=ncol(Input)
    modelF=TRUE
    }
  if(!CV){
    out=LARfit(Input,resp,q=size,modelFit = modelF)
    beta=out$beta
    error=out$RSS
    betaMat=out$betaMat
    errVector=out$errorVector
    yhat0=out$yhat0
    yhat=predictLS(X=Input,beta=beta)
  }
  if(CV & is.null(df)){
    source("D:/RProject/toolkits/ModelAssessment&Selection.R")
    DoF=seq(0,size)
    out=crossValidate(Input,resp,type=type,scaleInputs=cvScale,scaleType="Unit",complexity=DoF,NbreCV = 10)
    errVector=out$errorVector
    stdDev=out$stdDevVector
    error=min(errVector)
    bestSize=which.min(error)
    if(cvScale){Input=ScaleData(Input,resp,type = "Unit")$train}
    out00=LARfit(Input,resp,bestSize)
    beta=out00$beta
    yhat=predictLS(X=Input,beta=beta)
  }
  varList=which(beta[-1] !=0)

  return(list(yhat=yhat,beta=beta,RSS=error,betaMat=betaMat,varList=varList,
              errorVector=errVector,stdDev=stdDev,DoF=DoF))
}

#LAR fit
LARfit<-function(Input,resp,q,modelFit=FALSE,type="LAR"){
  resp0=resp-mean(resp)
  p=ncol(Input)
  N=nrow(Input)
  betaMat=matrix(0,ncol=q+1,nrow=p+1)
  betaMat[1,]=mean(resp)
  errorVector=rep(0,length=q+1)
  yhatMat=matrix(0,ncol=q+1,nrow=N)
  T_beta=c()
  S_beta=c()
  #q=0
  errorVector[1]=mean((resp-mean(resp))^2)
  yhatMat[,1]=mean(resp)
  
  if(q>0){
    #iteration 0
    yhat=0
    res=resp0-yhat

    for(i in 1:q){
      corrvector=t(Input)%*%res#covariates corrlation with residual
      if(i==1){#init C & active set
        C=max(abs(corrvector))
        Aj=which(abs(corrvector)==C)
      }

      sj=sign(corrvector[Aj])
      XA=as.matrix(Input[,Aj])%*%diag(sj)
      vv=c(10,20,30,40)
      if(ncol(XA)==1){GA=1}
      if(ncol(XA)>1){GA=solve(t(XA)%*%XA)}
      AA=1/sqrt(sum(GA))
      if(ncol(XA)==1){wA=1;uA=XA}
      if(ncol(XA)>1){
        wA=AA*as.matrix(rowSums(GA))
        uA=XA%*%wA #equi-angular vector
      }
      a=t(Input)%*%uA
      if(i<p){
        Ac=which(! seq(1:p) %in% Aj)
        gammaL=apply(as.matrix(Ac),1,function(x){
          tt=c((C-corrvector[x])/(AA-a[x]),(C+corrvector[x])/(AA+a[x]))
          vv=min(tt[tt>0])
          return(vv)
        })
        gamma=min(gammaL)
      }
      if(i==p){gamma=C/AA}
      betaMat[Aj+1,i+1]=betaMat[Aj+1,i]+sj*gamma*wA
      yhat=yhat+gamma*uA
      res=resp0-yhat
      T_beta=c(T_beta,sum(abs(betaMat[,i+1])))
      S_beta=c(S_beta,sum(res^2))
      yhatMat[,i+1]=yhat
      errorVector[i+1]=mean(res^2)
      if(i<q){#update C & active set for next step
        j0=Ac[which.min(gammaL)]
        Aj=c(Aj,j0)
        C=C-gamma*AA
      }
    }
  }
  if(modelFit){regSize=which.min(errorVector)}
  if(!modelFit){regSize=q+1}
  
  beta=betaMat[,regSize]
  yhat0=yhatMat[,regSize]
  RSS=errorVector[regSize]
  
  return(list(beta=beta,RSS=RSS,yhat0=yhat0,errorVector=errorVector,betaMat=betaMat,T_beta=T_beta,S_beta=S_beta))
}


#predictions for any linear model
predictLS<-function(LSModel=NULL,X,beta=NULL){
  if(!is.null(LSModel)){beta=LSModel$beta}
  yhat=as.matrix(cbind(1,X))%*%beta
  return(yhat)
}
