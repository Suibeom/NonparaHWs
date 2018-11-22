#This file contains functions relating to the log-likelihood
#of the normal exponential family, with a linear model for each
#of the two natural parameters.

lik=function(beta1,beta2,y,X1,X2){
  #y consists of n observations
  #X1 is an n-by-p1 matrix, with p1 covariates for each subject
  #X2 is an n-by-p2 matrix, with p2 covariates for each subject
  #beta1 is a vector of length p1 OR a m-by-p1 matrix
  #beta2 is a vector of length p2 OR a m-by-p2 matrix
  if (is.matrix(beta1)){
    beta1=t(beta1)
    beta2=t(beta2)
  }
  L=t(y)%*%X1%*%beta1 - t(y^2)%*%X2%*%beta2/2
  if (is.matrix(beta1)){
    L=L - 0.5*colSums((X1%*%beta1)^2/(X2%*%beta2))+0.5*colSums(log(X2%*%beta2))
    L=t(L)
    colnames(L)="Log-likelihood"
  } else {
    L=L - 0.5*sum((X1%*%beta1)^2/(X2%*%beta2))+0.5*sum(log(X2%*%beta2))
    L=L[1,1]
  }
  return(L-0.5*log(2*pi)*dim(X1)[1]) #returns a vector with m values of the likelihood
}

Dlik=function(beta1,beta2,y,X1,X2){
  #y consists of n observations
  #X1 is an n-by-p1 matrix, with p1 covariates for each subject
  #X2 is an n-by-p2 matrix, with p2 covariates for each subject
  #beta1 is a vector of length p1 OR a m-by-p1 matrix
  #beta2 is a vector of length p2 OR a m-by-p2 matrix
  
  if (is.matrix(beta1)){
    m=dim(beta1)[1]
    L=t(t(rep(1,m)))%*%c(t(y)%*%X1, - t(y^2)%*%X2/2)
    L=L - cbind(t((X1%*%t(beta1))/(X2%*%t(beta2)))%*%X1, -0.5*t((X1%*%t(beta1))^2/(X2%*%t(beta2))^2)%*%X2 -0.5*t(1/(X2%*%t(beta2)))%*%X2)   
  } else {
    L=c(t(y)%*%X1, - t(y^2)%*%X2/2)
    L=L - c(t((X1%*%beta1)/(X2%*%beta2))%*%X1, -0.5*t((X1%*%beta1)^2/(X2%*%beta2)^2)%*%X2 -0.5*t(1/(X2%*%beta2))%*%X2)   
  }
  if (is.matrix(beta1)) {
    colnames(L)=c(colnames(beta1),colnames(beta2))
    rownames(L)=rownames(beta1)} else {
      names(L)=c(names(beta1),names(beta2))
    }
  return(L) #returns m-by-(p1+p2) matrix with the m gradient vectors
}

D2lik=function(beta1,beta2,y,X1,X2){
  #y consists of n observations
  #X1 is an n-by-p1 matrix, with p1 covariates for each subject
  #X2 is an n-by-p2 matrix, with p2 covariates for each subject
  #beta1 is a vector of length p1 (only one parameter!)
  #beta2 is a vector of length p2 (only one parameter!)
  
  p1=length(beta1)
  p2=length(beta2)
  eta1=X1%*%beta1
  eta2=X2%*%beta2
  L=matrix(0,p1+p2,p1+p2)
  for (m in 1:(p1+p2)){
    for(n in 1:(p1+p2)){
      if (m<=p1 & n<=p1){
        L[m,n]=-sum(X1[,m]*X1[,n]/eta2)
      }
      if (m<=p1 & n>p1){
        L[m,n]=sum(eta1*X1[,m]*X2[,n-p1]/eta2^2)
      }
      if (m>p1 & n<=p1){
        L[m,n]=sum(eta1*X2[,m-p1]*X1[,n]/eta2^2)
      }
      if (m>p1 & n>p1){
        L[m,n]=-sum(X2[,m-p1]*X2[,n-p1]*(eta1^2/eta2^3 + 0.5/eta2^2))
      }
    }
  }
  #colnames(L)=c(names(beta1),names(beta2))
  #rownames(L)=c(names(beta1),names(beta2))
  
  return(L) #returns (p1+p2)-by-(p1+p2) second derivative matrix
}


