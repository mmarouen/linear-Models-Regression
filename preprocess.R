#PCA transformation
PCA<-function(Input){
  Input=as.matrix(Input)
  decomp=eigen(t(Input)%*%Input,symmetric = TRUE)
  D=diag(decomp$values)
  V=decomp$vectors
  Z=Input%*%V
  varProp=decomp$values/(sum(decomp$values));cumVarProp=cumsum(varProp)
  return(list(Z=Z,V=V,D=D,lamda=decomp$d,varianceProportion=varProp,cumulativeEnergy=cumVarProp))
}

ScaleData <- function(Input,Response=NULL,InputTest=NULL,ResponseTest=NULL,type=NULL,symm=TRUE,shuffle=FALSE){
  N=nrow(Input)
  means=c()
  stds=c()
  if(! is.null(Response) & shuffle){
    inds=sample(N,N)
    Input=Input[inds,]#avoid DB prior filtering
    Response=Response[inds]
  }
  if((! is.null(InputTest))&(! is.null(ResponseTest))& shuffle){
    N1=nrow(InputTest)
    inds1=sample(N1,N1)
    InputTest=InputTest[inds1,]
    ResponseTest=ResponseTest[inds1]
  }
  if(!is.null(type)){
    #extract only quantitative inputs
    if(ncol(Input)>1){
      cols=!apply(Input,2,is.factor)
      tmp=Input[,cols]
    }
    if(ncol(Input)==1){
      tmp=as.matrix(Input)
      cols=1
    }
    if(type=="Unit"){
      tmp=scale(tmp,TRUE,TRUE)
      means=attr(tmp,"scaled:center")
      stds=attr(tmp,"scaled:scale")
      tmp=apply(tmp,2,function(x) x/sqrt(sum(x^2)))
      Input[,cols]=tmp
      if(! is.null(InputTest)){
        tmp1=as.matrix(InputTest[,cols])
        tmp1=t(apply(tmp1,1,'-', means ))
        tmp1=t(apply(tmp1,1,'/', stds ))
        tmp1=apply(tmp1,2,function(x) x/sqrt(sum(x^2)))
        InputTest[,cols]=as.matrix(tmp1)
      }
    }
    if(type=="Standardize"){#studentized residuals
      tmp=scale(tmp,TRUE,TRUE)
      means=attr(tmp,"scaled:center")
      stds=attr(tmp,"scaled:scale")
      Input[,cols]=tmp
      if(! is.null(InputTest)){
        tmp1=as.matrix(InputTest[,cols])
        tmp1=t(apply(tmp1,1,'-', means ))
        tmp1=t(apply(tmp1,1,'/', stds ))
        InputTest[,cols]=as.matrix(tmp1)
      }
    }
    if(type=="Normalize"){#covariates in range 0,1 or -1,1
      Input[,cols]=apply(tmp,2,function(x) (x-min(x))/(max(x)-min(x)))
      if(! is.null(InputTest)){
        tmp1=InputTest[,cols]
        InputTest[,cols]=apply(tmp1,2,function(x) (x-min(x))/(max(x)-min(x)))
      }
      if(symm){
        Input=2*Input-1
        if(! is.null(InputTest)){InputTest=2*InputTest-1}
      }
    }
  }
  
  return(list(train=Input,respTrain=Response,test=InputTest,respTest=ResponseTest,paramList=list(means=means,stds=stds)))
}

