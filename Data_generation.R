#Generating data


#packages
library(parallel)
library(foreign)
library(Rlab)
library(stats4)
library(mcmc)
library(MCMCpack)
library(R2jags)
library(plyr)
library(dplyr)
library(R2WinBUGS) 
library("fdrtool")   
library(AER)
library(akima)

#simulation settings
sim=seed=1                 #Use the same seed for each scenario
Niter=10000                #number of mcmc samples
n_patients<- 100           #Number of patients in simulated data set
H=n_hist_trials <- 2       #Number of historical trials
var_study_effect <- 0
lambdac <- 4               #average rate in the control arm
baseline_linear <- log(4  )#Log baseline odds, log(0.72/0.28) from binomoal
r=0.4                      #These value depend on scenarios
r1=Inf  
r2= 0.8
intervention_effect<- -0.6 #Treatment effect 


#Simulate the data set
trialnr<-rep(0,(n_hist_trials+2)*n_patients)
response<- rep(0,(n_hist_trials+2)*n_patients)
intervention<- rep(0,(n_hist_trials+2)*n_patients)

#Historical trials
for (i in 1:(n_hist_trials)){
  #set.seed(2)
  trial_effect<-rnorm(n=1,mean=0,sd=sqrt(var_study_effect)) + baseline_linear
  trial_effect_indiv<- trial_effect
  for (j in 1:n_patients){    
    trialnr[(i-1)*n_patients+j]<-i
    if(i==1){
      response[(i-1)*n_patients+j]<- rnbinom(1,size=r1,mu=exp(trial_effect_indiv))}
    if(i==2){
      response[(i-1)*n_patients+j]<- rnbinom(1,size=r2,mu=exp(trial_effect_indiv))}  
  }  
}
#Current trial
i<- n_hist_trials+1
random_trial_effect_current_trial<-rnorm(n=1,mean=0,sd=sqrt(var_study_effect))
for (j in 1:n_patients){
  #Control group of current trial
  trial_effect<-random_trial_effect_current_trial+baseline_linear
  trialnr[(i-1)*n_patients+j]<-n_hist_trials+1
  trial_effect_indiv<- trial_effect
  response[(i-1)*n_patients+j]<-rnbinom(1,size=r,mu=exp(trial_effect_indiv))
} 
i<-n_hist_trials+2
for (j in 1:n_patients){
  #Intervention group of current trial
  trial_effect<-random_trial_effect_current_trial+baseline_linear+intervention_effect
  intervention[(i-1)*n_patients+j]<-1
  trialnr[(i-1)*n_patients+j]<- n_hist_trials+1
  trial_effect_indiv<- trial_effect
  response[(i-1)*n_patients+j]<-rnbinom(1,size=r,mu=exp(trial_effect_indiv))
}  
dataset2<-data.frame(trialnr=trialnr,response=response,Intervention=intervention)

head(dataset2)
Y0=list(H)
Y0[[1]]<- dataset2[dataset2$trialnr==1,]$response; length(Y0[[1]])
Y0[[2]]<- dataset2[dataset2$trialnr==2,]$response; length(Y0[[2]])
Y00=c(Y0[[1]],Y0[[2]])
N0=length(Y00)
aaa=dataset2[dataset2$trialnr==3,]
Yc=aaa[aaa$Intervention==0,]$response
Yt=aaa[aaa$Intervention==1,]$response
nc=length(Yc)
nt=length(Yt)
Y=c(Yc,Yt)
Yc_pooled=c(Y0[[1]],Y0[[2]],Yc)
Nc= length(Yc_pooled)
Y_pooled=c(Y0[[1]],Y0[[2]],Y)
N=length(Y_pooled)
Nct=nc+nt

X0=list(n_hist_trials)
X0[[1]]=cbind(1,dataset2[dataset2$trialnr==1,]$Intervention)
X0[[2]]=cbind(1,dataset2[dataset2$trialnr==2,]$Intervention)
Xc=cbind(1,aaa[aaa$Intervention==0,]$Intervention)
Xt=cbind(1,aaa[aaa$Intervention==1,]$Intervention)
X=rbind(Xc,Xt)
X_pooled=rbind(X0[[1]],X0[[2]],X)
dim=ncol(X0[[1]]) 
Trial<-as.numeric(dataset2$trialnr[1:N])

n0=rep(0,H)
slyf0=rep(0,H)
s0=rep(0,H)
for(i in 1:(H)){
  s0[i]= sum(Y0[[i]])
  slyf0[i]= sum(log(factorial(Y0[[i]])))
  n0[i]=length(Y0[[i]])
}
sc= sum(Yc)
slyfc= sum(log(factorial(Yc)))
st= sum(Yt)
slyft= sum(log(factorial(Yt)))
HistStudy<-as.numeric(dataset2$trialnr[1:N0])
sigma_data= 2  # not sd(Y), It is defined as sqrt(k+1/mean(Y)) in Gstieger et al 
precision_tau<-1/((sigma_data/2)^2);precision_tau # precision of between-trial-heterogeniety
precision2<-1/((sigma_data)^2);precision2 # precision of the robust component





