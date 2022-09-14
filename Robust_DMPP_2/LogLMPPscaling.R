

#Only for calculating the scaling constant, so current trial is not included in calculation!
LikelihoodMPPscaling<- function(parameters,weight1,weight2,X0,Y0){
  dim=ncol(X0[[1]])
  beta <- parameters[1:dim]
  k=parameters[dim+1] 
  r<- 1/k
  if(weight1>=0 & weight1<=1&weight2>=0 & weight2<=1 &r>0){ 
    lambda01 <- as.vector(exp(X0[[1]]%*% beta))
    lambda02 <- as.vector(exp(X0[[2]]%*% beta))
    lik.01<-   weight1*sum(lgamma(Y0[[1]]+r) - lfactorial(Y0[[1]]) - lgamma(r) +r * log(r/(r+lambda01)) + Y0[[1]] * log(lambda01/(r+lambda01)))
    lik.02<-   weight2*sum(lgamma(Y0[[2]]+r) - lfactorial(Y0[[2]]) - lgamma(r) +r * log(r/(r+lambda02)) + Y0[[2]] * log(lambda02/(r+lambda02)))
    Prior0<-dnorm(beta[1],0,5,log=TRUE)+dnorm(beta[2],0,2,log=TRUE)
    Prior0<-Prior0+dexp(k,0.4,log=TRUE)
  }else{
    lik.01<--10^200
    lik.02<--10^200
    Prior0 <- -10^200
  }
  LikelihoodMPP_Ind<- Prior0+lik.01+lik.02
}

#Log-likelihood function (excluding prior) for MPP.
#Only for calculating the scaling constant, so current trial is not included in calculation.
LikelihoodMPPwithoutprior <- function(parameters,X0,Y0){
  dim=ncol(X0[[1]])
  beta <- parameters[1:dim]
  k=parameters[dim+1] 
  r<- 1/k
  Y01=as.vector(Y0[[1]])
  Y02=as.vector(Y0[[2]])
  lambda01 <- as.vector(exp(X0[[1]]%*% beta))
  lambda02 <- as.vector(exp(X0[[2]]%*% beta))
  lik.01<-    sum(lgamma(Y0[[1]]+r) - lfactorial(Y0[[1]]) - lgamma(r) +r * log(r/(r+lambda01)) + Y0[[1]] * log(lambda01/(r+lambda01)))
  lik.02<-    sum(lgamma(Y0[[2]]+r) - lfactorial(Y0[[2]]) - lgamma(r) +r * log(r/(r+lambda02)) + Y0[[2]] * log(lambda02/(r+lambda02)))
  LikelihoodMPP_Ind<- list(Loglik1=lik.01,Loglik2=lik.02)
}

