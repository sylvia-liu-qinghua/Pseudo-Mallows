rm(list=ls())
require(BayesMallows)
n<-200  
N<-200
rho0<-1:n
runs<-1:2
alpha0s<-c(3)
lambdas<-c(5)
load("./shared/Cdfootrule.RData")
source("./shared/pseudoSamplingFuncs.R")
if(n>50){
  fitvec = estimate_partition_function(alpha_vector = seq(0.01,10,0.2), n_items = 50,metric = "footrule", nmc = 2000,degree=10)
}
for(alpha0 in alpha0s){
  for(lambda in lambdas){
    for(run in runs){
      origin_data<-sample_mallows(rho0 = rho0,alpha0 = alpha0, n_samples = N)
      n<-dim(origin_data)[2]
      N<-dim(origin_data)[1]
      clicking_data<-createClickData(origin_data,lambda)
      save(origin_data,clicking_data,file=paste("./datasets/N",N,'n',n,'alpha',alpha0,'lambda',lambda, 'run',run,'.RData',sep=""))
      
    }
  }
}
