

#Log-likelihood function for Robust_DMPP_2
LikelihoodMPP_Robust_DMPP_2 <- function(parameters,H,X0,X,Y0,Y,weightTemp,ScalingConstantTemp){
  dim=ncol(X0[[1]])
  beta <- c(parameters[1],parameters[2])
  weight=parameters[(dim+1):(dim+H)] 
  k=parameters[dim+H+1] 
  a=parameters[dim+H+2] 
  b=parameters[dim+H+3]
  r<- 1/k
  mu=a/(a+b)
  sigmasq=(mu*(1-mu))/(a+b+1)
  Y01=as.vector(Y0[[1]])
  Y02=as.vector(Y0[[2]])
  if(any(weight>=0 & weight<=1) &r>0& mu>0 & mu <1 & sigmasq>0 & sigmasq<1 & a>0 & b>0){ 
    lambda01 <- as.vector(exp(X0[[1]]%*% beta))
    lambda02 <- as.vector(exp(X0[[2]]%*% beta))
    lik.01<-   weight[1]*sum(lgamma(Y01+r) - lfactorial(Y01) - lgamma(r) +r * log(r/(r+lambda01)) + Y01 * log(lambda01/(r+lambda01)))
    lik.02<-   weight[2]*sum(lgamma(Y02+r) - lfactorial(Y02) - lgamma(r) +r * log(r/(r+lambda02)) + Y02 * log(lambda02/(r+lambda02)))
    lambda <- as.vector(exp(X%*% beta))
    lik<-      sum(lgamma(Y+r) - lfactorial(Y) - lgamma(r) +r * log(r/(r+lambda)) + Y * log(lambda/(r+lambda)))
    Prior0<-dnorm(beta[1],0,5,log=TRUE)+dnorm(beta[2],0,2,log=TRUE)
    Prior0<- Prior0+log(dunif(mu,0,1))
    Prior0<-Prior0+log(dgamma(1/sigmasq,0.01,0.01)/(sigmasq^2))
    Prior0<-Prior0+log (0.9*prod(dbeta(weight,a,b))+0.1*prod(dhalfnorm(weight,sd2theta(sqrt(sigmasq/6.25)))))
    Prior0<-Prior0+dexp(k,0.4,log=TRUE)
    ScalingConstantTemp <- bilinear(weightTemp, weightTemp, ScalingConstantTemp, weight[1], weight[2])$z  
  }else{
    lik.01<--10^200
    lik.02<--10^200
    lik<- -10^200
    Prior0 <- -10^200
    ScalingConstantTemp<- 0
  }
  LikelihoodMPP_Robust_DMPP_2<- Prior0+lik.01+lik.02+lik-ScalingConstantTemp
}


