#compute expected prediction error using cross validation N-folds
#supported: PCR, ridge, LS, PLS, Lasso, LAR, IFSR, Neural Nets, regularized DA, reduced rank DA
crossValidate <- function(X,resp,type,complexity=NULL,NbreCV=10,
                          scaleInputs=FALSE,scaleType=NULL,symm=FALSE){
  shuffleData=ScaleData(X,resp,shuffle = TRUE)
  X=shuffleData$train;resp=shuffleData$respTrain
  X=as.matrix(X)
  N=nrow(X)
  p=ncol(X)
  ntest=round(N/NbreCV)
  ntrain=N-ntest
  errorVector=rep(0,length = NbreCV)
  sdVal=c()
  errorMat=matrix(0,ncol=NbreCV,nrow = length(complexity))
  sdVector=rep(0,length=length(complexity))
  yHatTest=c()
  yHatTestMat=c()
  betaMat=c()
  for (i in 1:NbreCV){
    #subset CV data
    testInds=(i-1)*ntest+1:ntest
    testInds=intersect(testInds,1:N)
    trainInds=setdiff(1:N,testInds)
    train_CV=as.matrix(X[trainInds,])
    resptrain=resp[trainInds]
    test_CV=as.matrix(X[testInds,])
    resptest=resp[testInds]
    #preprocessing
    if(scaleInputs){
      source("D:/RProject/toolkits/Preprocess&FS.R")
      cvScale=ScaleData(train_CV,Response = resptrain,InputTest = test_CV,type = scaleType,symm = symm,shuffle = TRUE)
      train_CV=cvScale$train;test_CV=cvScale$test;resptrain=cvScale$respTrain
    }
    #compute yHat for test data
    if(type == "LS"){
      source("D:/RProject/toolkits/Regression.R")
      out=LS(train_CV,resptrain,intercept = complexity)
      yHatTest=predictLS(out,test_CV)
    }
    if(type == "PCR"){
      source("D:/RProject/toolkits/Regression.R")
      q=ncol(train_CV)
      out=PCRfit(train_CV,resptrain,q = q)
      betaMat=out$betaMat
      yhatTestMat=apply(betaMat,2,function(x) predictLS(X = test_CV,beta = x))
    }
    if(type == "Ridge"){
      source("D:/RProject/toolkits/Regression.R")
      out=ridgefit(train_CV,resptrain,lambda = complexity)
      yHatTest=predictLS(X = test_CV,LSModel = out)
    }
    if(type=="PLS"){
      source("D:/RProject/toolkits/Regression.R")
      q=complexity[length(complexity)]
      out=PLSfit(train_CV,resptrain,q)
      betaMat=out$betaMat
      yhatTestMat=apply(betaMat,2,function(x) as.matrix(cbind(1,test_CV))%*%x)
    }
    if(type == "RRDA"){
      out=RRDAfit(train_CV,resptrain,test_CV,modelF=FALSE,sub=complexity)
      yHatTest=out$yHatTest
    }
    if(type=="RegDA"){
      alpha=complexity[1]
      gamma=complexity[2]
      out=RegDAfit(train_CV,resptrain,InputTest=test_CV,alpha,gamma)
      yhatTestMat=out$yHatTest
    }
    if(type=="RF"){
      train_CV=as.data.frame(train_CV)
      library(randomForest)
      resptrain=factor(resptrain)
      levels(resptrain)=levels(resp)
      out11=randomForest(resptrain~.,data=cbind(resptrain,train_CV))
      yHatTest=as.numeric(as.character(predict(out11,test_CV)))
    }
    if(type=="NN"){
      out11=neuralNet(train_CV,resptrain,DW=c(12),type="Classification",loss="Deviance",outputFunc = "Softmax",activationFunc = "tanh",rate=complexity[[1]],weightDecay = complexity[[2]],lambda = complexity[[3]])
      yHatTest=predictNN(out11,test_CV)
    }
    if(type=="Subset"){
      source("D:/RProject/toolkits/Regression.R")
      out=subsetfit(train_CV,resptrain,q = complexity)
      yHatTest=predictLS(X = test_CV,LSModel = out)
    }
    if(type %in%c("Lasso","LAR","IFSR")){
      source("D:/RProject/toolkits/Regression.R")
      q=ncol(train_CV)
      out=LARfit(train_CV,resptrain,q = q,modelFit = FALSE,type = type)
      betaMat=out$betaMat
      yhatTestMat=apply(betaMat,2,function(x) predictLS(X = test_CV,beta = x))
    }
    
    #prediction error compute
    if(type %in% c("LS","Ridge","Subset")){
      errorVector[i]=mean((yHatTest-resptest)^2)}
    if(type %in% c("RRDA","RegDA","NN","RF")){
      errorVector[i]=mean(as.numeric(as.character(yHatTest)) != resptest)}
    if(type %in% c("PLS","IFSR","Lasso","LAR","PCR")){
      errorMat[,i]=apply(yhatTestMat,2,function(x) mean((x-resptest)^2))
    }
  }
  if(type %in% c("PLS","IFSR","Lasso","LAR","RRDA","PCR")){
    errorVector=rowMeans(errorMat)
    }
  sdVal=sd(errorVector)/sqrt(NbreCV)
  sdVector=apply(errorMat,1,sd)*(1/sqrt(NbreCV))
  return(list(errorVector=errorVector,stdDev=sdVal,stdDevVector=sdVector,betaMat=betaMat))
}

#CV plotting with standard deviation for error values
CVplot<-function(errVector,stdVector,DoF,stepSize=1){
  epsilon = 0.05
  L=length(errVector)#must be same length for DoF and stdVector
  pins=seq(1,L,by=stepSize)
  L1=length(pins)
  errVec=errVector[pins]
  stdVec=stdVector[pins]
  Df=DoF[pins]
  up=errVec+stdVec
  low=errVec-stdVec
  plot.new()
  plot(Df,errVec,type = 'b',pch=19,ylim=c(min(low),max(up)),xlab="Degrees of freedem",ylab = "CV error",col="orange")
  for(i in 1:L1) {
    segments(Df[i],low[i],Df[i],up[i],col="blue")
    segments(Df[i]-epsilon, up[i] ,Df[i]+epsilon, up[i],col="blue")
    segments(Df[i]-epsilon, low[i] ,Df[i]+epsilon, low[i],col="blue")
  }
}

#bootstrap used for model selection (predict the EPE) using 03 methods
# 1.1 Standard: error is computed on the whole training set (out-sample training error)
# 1.2 ".632": a compromise between in-sample training error & prediction error
# 1.3 LOO (leave one out): compute mean of generalization errors
#supported models: PCR, ridge, LS, Lasso, LAR, IFSR, Neural Nets,reduced rank DA
#not yet supported: PLS, Regularized discriminant analysis
bootstrap<-function(Input,resp,type,B=100,variation="LOO",scaleInputs=FALSE,
                    scaleType=NULL,symm=FALSE,intercept=FALSE,loss="RSS"){
  source("D:/RProject/toolkits/Regression.R")
  source("D:/RProject/toolkits/Preprocess&FS.R")
  
  Input=as.matrix(Input)
  N=nrow(Input)
  errorMat=matrix(-1,ncol=B,nrow=N)
  errBoot=rep(0,length=B)
  errLOO=c()
  #init datasets
  X=Input#used to sample from
  Xtest=X#used as testing set
  yHatTestMat=c()
  yHatTest=c()
  err1=c()
  errboot_enhanced=c()
  for(b in 1:B){
    inds_b=sample(N,N,replace=TRUE)
    Xb=X[inds_b,]
    yb=resp[inds_b]
    #preprocess
    if(scaleInputs){
      bootScale=ScaleData(Xb,Response = yb,InputTest = Xtest,type = scaleType,symm = symm)
      Xb=bootScale$train;Xtest=bootScale$test;yb=bootScale$respTrain
    }
    
    #model
    if(type=="Customize"){
      Hb=bs(Xb,knots = c(0.75,1.5,2.25),df = 7,intercept = TRUE)
      Xtestb=bs(Xtest,knots = c(0.75,1.5,2.25),df = 7,intercept = TRUE)
      out=LS(Hb,yb,inputTest=Xtestb,scaleIn=scaleInputs)
      yHatTest=out$yhatTest
    }
    if(type == "LS"){
      out=LR(Xb,yb)
      beta=out$beta
      yHatTest=as.matrix(cbind(1,Xtest))%*%beta
    }
    if(type == "PCR"){
      out=PCR(Xb,yb,InputTest=Xtest,size=complexity)
      betaHat=out$beta
      yHatTest=out$yhat
    }
    if(type=="RF"){
      Xb=as.data.frame(Xb)
      library(randomForest)
      yb=factor(yb)
      out11=randomForest(yb~.,data=cbind(yb,Xb))
      yHatTest=as.numeric(as.character(predict(out11,Xtest)))
    }
    if(type=="LDA"){
      source("D:/RProject/toolkits/Classification.R")
      Xb=as.data.frame(Xb)
      Xb$AGE=log(Xb$AGE+0.1);Xb$ARPers=log(Xb$ARPers+0.1)
      Xt=as.data.frame(Xtest)
      Xt$AGE=log(Xt$AGE+0.1);Xt$ARPers=log(Xt$ARPers+0.1)
      # source("D:/RProject/toolkits/Preprocess&FS.R")
      # scaled=ScaleData(Input = Xb,Response = yb,InputTest = Xt,shuffle = TRUE,type = "Standardize")
      # Xb=scaled$train
      # yb=scaled$respTrain
      # Xt=scaled$test
      fit1=LinDA(Xb,yb)
      yHatTest=predictLinDA(fit1,Xt)
    }
    if(type=="GBM"){
      library(xgboost)
      xgbmfit=xgboost(data = as.matrix(Xb),label = yb,nrounds = 5000,verbose = 0)
      yy=predict(xgbmfit,as.matrix(Xb))
      yHatTest=round(yy-1)
    }
    if(type=="LR"){
      library(nnet)
      Xb=as.data.frame(Xb)
      Xb$AGE=log(Xb$AGE+0.1);Xb$ARPers=log(Xb$ARPers+0.1)
      Xt=as.data.frame(Xtest)
      Xt$AGE=log(Xt$AGE+0.1);Xt$ARPers=log(Xt$ARPers+0.1)
      source("D:/RProject/toolkits/Preprocess&FS.R")
      scaled=ScaleData(Input = Xb,Response = yb,InputTest = Xt,shuffle = TRUE,type = "Standardize")
      Xb=scaled$train
      yb=scaled$respTrain
      Xt=scaled$test
      fit1=multinom(yb~.,data=cbind(Xb,yb),trace=FALSE)
      yHatTest=as.numeric(as.character(predict(fit1,newdata = Xt,type = "class")))
    }
    if(type == "Ridge"){
      out=ridge(Xb,yb,Xtest,complexity)
      yHatTest=out$yhat
      betaHat=out$beta
    }
    if(type=="KNN"){
      library(class)
      yHatTest=knn(Xb,Xtest,cl = yb,k = 10)
    }
    if(type == "RRDA"){
      out=RRDA(Xb,yb,Xtest)
      yHatTest=out$yHatTest
    }
    if(type=="NN"){
      out11=neuralNet(Xb,yb,DW=c(5),type="Classification",loss="Deviance",outputFunc = "Softmax",activationFunc = "tanh",rate=complexity[[2]],weightDecay = complexity[[3]],lambda = complexity[[4]])
      yHatTest=predictNN(out11,Xtest)
    }
    if(type %in%c("Lasso","LAR","IFSR")){
      out=LARS(Xb,yb,type,t=complexity)
      if(is.null(out$covariates)){yHatTest=mean(yb)}
      if(!is.null(out$covariates)){
        Xtest=apply(Xtest,2,function(x) x/sqrt(sum(x^2)))
        yHatTest=Xtest%*%out$beta+mean(yb)
      }
    }
    
    eps=rnorm(N,0,1)
    yHatTestMat=cbind(yHatTestMat,yHatTest)
    
    #compute error by bootstrap sample b
    if(variation == "boot"){
      if(loss=="RSS"){
        errBoot[b]=mean((resp-yHatTest)^2)
      }
      if(loss=="Classification"){
        errBoot[b]=mean(resp != yHatTest)
      }
    }
    if(variation != "boot"){
      idx=!seq(1,N)%in%inds_b
      if(loss=="RSS"){
        errorMat[idx,b]=(resp[idx]-yHatTest[idx])^2
        errBoot[b]=mean(errorMat[idx,b])
      }
      if(loss=="Classification"){
        errorMat[idx,b]=resp[idx] != yHatTest[idx]
        errBoot[b]=mean(resp[idx] != yHatTest[idx])
        #print(errBoot[b])
      }
    }
  }
  if(variation == ".632"){#we need training error
    Xtrain=X
    if(type == "LS"){
      betaHat=solve(t(Xtrain)%*%Xtrain)%*%t(Xtrain)%*%(resp-mean(resp))
      yHatTrain=Xtest%*%betaHat+mean(yb)
    }
    if(type == "PCR"){
      out=PCR(Xtrain,resp,InputTest=Xtrain,size=complexity)
      betaHat=out$beta
      yHatTrain=out$yhat
    }
    if(type == "Ridge"){
      out=ridge(Xtrain,resp,Xtrain,complexity)
      yHatTrain=out$yhat
      betaHat=out$beta
    }
    if(type == "RRDA"){
      out=RRDA(Xtrain,resp,Xtrain)
      yHatTrain=out$yHatTest
    }
    if(type=="RF"){
      Xtrain=as.data.frame(Xtrain)
      resp=factor(resp)
      out11=randomForest(resp~.,data=cbind(resp,Xtrain))
      yHatTrain=as.numeric(as.character(predict(out11,Xtrain)))
    }
    if(type=="GBM"){
      library(xgboost)
      xgbmfit=xgboost(data = as.matrix(Xtrain),label = resp,nrounds = 5000,verbose = 0)
      yy=predict(xgbmfit,as.matrix(Xtrain))
      yHatTrain=round(yy-1)
    }
    if(type=="LDA"){
      source("D:/RProject/toolkits/Classification.R")
      Xtrain=as.data.frame(Xtrain)
      Xtrain$AGE=log(Xtrain$AGE+0.1);Xtrain$ARPers=log(Xtrain$ARPers+0.1)
      # scaled=ScaleData(Input = Xtrain,Response = resp,shuffle = TRUE,type = "Standardize")
      # Xtrain=scaled$train
      # yTrain=scaled$respTrain
      fit1=LinDA(Xtrain,resp)
      yHatTrain=predictLinDA(fit1,Xtrain)
    }
    if(type=="LR"){
      library(nnet)
      source("D:/RProject/toolkits/Preprocess&FS.R")
      Xtrain=as.data.frame(Xtrain)
      Xtrain$AGE=log(Xtrain$AGE+0.1);Xtrain$ARPers=log(Xtrain$ARPers+0.1)
      scaled=ScaleData(Input = Xtrain,Response = resp,shuffle = TRUE,type = "Standardize")
      Xtrain=scaled$train
      yTrain=scaled$respTrain
      fit1=multinom(yTrain~.,data=cbind(Xtrain,yTrain),trace=FALSE)
      yHatTrain=as.numeric(as.character(predict(fit1,newdata = Xtrain,type = "class")))
    }
    if(type=="KNN"){
      library(class)
      yHatTrain=knn(Xtrain,Xtrain,cl = yb,k = 10)
    }
    if(type=="NN"){
      out11=neuralNet(Xtrain,resp,DW=c(5),type="Classification",loss="Deviance",outputFunc = "Softmax",activationFunc = "tanh",rate=complexity[[2]],weightDecay = complexity[[3]],lambda = complexity[[4]])
      yHatTrain=predictNN(out11,Xtrain)
    }
    if(type %in%c("Lasso","LAR","IFSR")){
      out=LARS(Xtrain,resp,type,t=complexity)
      if(is.null(out$covariates)){yHatTrain=mean(resp)}
      if(!is.null(out$covariates)){
        Xtrain=apply(Xtrain,2,function(x) x/sqrt(sum(x^2)))
        yHatTrain=Xtrain%*%out$beta+mean(resp)
      }
    }
  }
  
  #overall bootstrap error
  if(variation=="boot"){err1=mean(errBoot)}
  if(variation != "boot"){ 
    
    #### compute out of the bag error
    errLOO=apply(errorMat,1,function(x){
      C_minus_i= x!= -1
      err_minus_i=mean(x[C_minus_i])
      return(err_minus_i)
    })
    
    if(variation==".632"){
      
      if(loss=="RSS"){
        trainErr=mean((resp-yHatTrain)^2)
      }
      if(loss=="Classification"){
        trainErr=mean(resp != yHatTrain)
      }
      err1=0.368*trainErr+0.632*mean(errLOO)
      errboot_enhanced=0.632*errBoot+0.368*trainErr
    }
  }
  return(list(bootstrapError=err1,yHatTest=yHatTestMat,errBoot=errBoot,errboot1=errboot_enhanced))
}

