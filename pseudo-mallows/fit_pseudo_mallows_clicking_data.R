rm(list=ls())
require(RcppArmadillo)
require(Rcpp)
require(BayesMallows)
require(fields)
# setwd("~/Desktop/Oslo/Spring2020/Simulation_Partial/")

sourceCpp("./Cpp_R_Implementation/clickingAlgo.cpp")
source("./shared/allFunctions.R")
source("./shared/pseudoSamplingFuncs.R")
runs<-1
alpha0s<-c(3)
n<-200
N<-200
rho0<-1:n
lambdas<-c(5)
n_samples<-c(100,500,1000)
for(alpha0 in alpha0s){
  for(run in runs){
    for(lambda in lambdas){
      for(n_sample in n_samples){
        print(paste("alpha0 = ", alpha0))
        print(paste("run = ", run))
        load(paste("./datasets/N",N,'n',n,'alpha',alpha0,'lambda',lambda, 'run',run,'.RData',sep=""))
        centre_inferred<-n+1-rank(apply(clicking_data,2,sum),ties.method = "first")
        data_init<-t(apply(clicking_data,1,generateInit_oneuser_2,centre = centre_inferred))
        
        start = proc.time()
        results1<-PseudoForClicks(n_samples=n_sample,clicking_data=clicking_data, data_init = data_init, alpha =alpha0*2,sigma=0)
        durationCpp1= proc.time()-start
        
        K<-10
        userArray3D_cpp1<-results1$individuals[,,(n_samples/2+1):n_samples]
        recsCpp1<-vector()
        truthCpp<-vector()
        
        for(userj in 1:N){
          print(userj)
          heatMat_userjCpp1<-heatMap(t(userArray3D_cpp1[userj,,]),1:n)
          recsCpp1<-rbind(recsCpp1,recommendK(heatMat_userjCpp1,K, clicking_data[userj,],type="click"))
          truthCpp<-rbind(truthCpp, origin_data[userj,][origin_data[userj,]%in%((sum(clicking_data[userj,])+1):(sum(clicking_data[userj,])+K))])
        }
        
        correctOrWrongCpp1<-recsCpp1
        for(userj in 1:N){
          correctOrWrongCpp1[userj,]<-recsCpp1[userj,]%in%truthCpp[userj,]
        }
        save(n_sample,durationCpp1, file=paste("./results/durations/N",N,'n',n,'alpha',alpha0,'lambda',lambda, 'run',run,'nmc',n_sample,'.RData',sep=""))
        save(results1, recsCpp1, correctOrWrongCpp1,file=paste("./results/Recs/N",N,'n',n,'alpha',alpha0,'lambda',lambda, 'run',run,'nmc',n_sample, '.RData',sep=""))
      }  
    }
  }
}
    