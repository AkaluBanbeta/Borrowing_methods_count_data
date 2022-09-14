
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

