
niter=1000             #number of iteration in path sampling
gridsize=0.05          #step size of the powers from 0 and 1

gweight1 <- seq(0,1,gridsize)+0.5*gridsize
gweight1<-gweight1[1:(length(gweight1)-1)]
gridsize0=gridsize/4
gweight0=seq(0,0.05,gridsize0)+0.5*gridsize0
gweight=c(gweight0,gweight1[gweight1>max(gweight0)])

deriv_scaling_constant_alpha1<-matrix(NA,nrow=length(gweight),ncol=length(gweight)) 
deriv_scaling_constant_alpha2<-matrix(NA,nrow=length(gweight),ncol=length(gweight)) 
output_deriv_alpha1<-rep(0,niter/2)
output_deriv_alpha2<-rep(0,niter/2)
init1s=rep(1,dim+1)

for(index_weight2 in 1:length(gweight)){
  for(index_weight1 in 1:length(gweight)){
    weight1<-gweight[index_weight1]
    weight2<-gweight[index_weight2]
    v=sqrt(0.01/weight2);print(v) # step size of the metropolis hastings
    print(weight1)
    print(weight2)
    if (index_weight1==1){ 
      result_MPP_temp<-MCMCmetrop1R(LikelihoodMPPscaling,theta.init=init1s,burnin = 1000,mcmc =niter,
                                    seed=1,verbose=0,thin = 1 ,V=v*diag(length(init1s)),weight1=weight1,
                                    weight2=weight2,X0=X0,Y0=Y0)
    }
    if(index_weight1>1){
      result_MPP_temp<-MCMCmetrop1R(LikelihoodMPPscaling,theta.init=result_MPP_temp[nrow(result_MPP_temp),],
                                    seed=1,burnin = 0,mcmc = niter,verbose=0, thin = 1 ,V=1*diag(pmax(0.00005,
                                                                                                      diag(cov(result_MPP_temp)))), weight1=weight1,weight2=weight2,X0=X0,Y0=Y0)
    }
    for (ii in 1:(nrow(result_MPP_temp)/2)){ # what if this is 'for (ii in 1:(nrow(result_MPP_temp)/2)){'
      res<-LikelihoodMPPwithoutprior(result_MPP_temp[(ii+nrow(result_MPP_temp)/2),],X0=X0,Y0=Y0)
      output_deriv_alpha1[ii]<- res$Loglik1
      output_deriv_alpha2[ii]<- res$Loglik2
    } 
    deriv_scaling_constant_alpha1[index_weight1,index_weight2]<-mean(output_deriv_alpha1)
    deriv_scaling_constant_alpha2[index_weight1,index_weight2]<-mean(output_deriv_alpha2)
  }
}

rowCumSums <- function(x) {
  for(i in seq_len(dim(x)[1])) { x[i,] <- cumsum(x[i,]) }; x
}
colCumSums <- function(x) {
  for(i in seq_len(dim(x)[2])) { x[,i] <- cumsum(x[,i]) }; x
}

###########  calculation of partial deriv_scaling_constant for alpha1 and alpha2
weightTemp=c(0,gweight0+0.5*gridsize0,gweight[(length(gweight0)+1):length(gweight)]+0.5*gridsize) 
# actual weights in analysis, starting from 0
stepsize=diff(weightTemp) # the stepsize of the weight values
deriv_scaling_constant_alpha12=deriv_scaling_constant_alpha1*stepsize
deriv_scaling_constant_alpha22=t(t(deriv_scaling_constant_alpha2)*stepsize)


########### 2. calculation of scaling constant 
ScalingConstant<-0.5*rowCumSums(deriv_scaling_constant_alpha22)+ 
  0.5*rowCumSums(deriv_scaling_constant_alpha22)[c(2:length(gweight),length(gweight)),]+ 
  1.5*colCumSums(deriv_scaling_constant_alpha12)[,rep(1,length(gweight))]-
  0.5* colCumSums(deriv_scaling_constant_alpha12)[,rep(2,length(gweight))]
ScalingConstantTemp=rbind(0,cbind(0,ScalingConstant))

