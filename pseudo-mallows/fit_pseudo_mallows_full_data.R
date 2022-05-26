rm(list=ls())
require(RcppArmadillo)
require(Rcpp)
require(BayesMallows)
require(fields)
setwd("~/Desktop/Oslo/Spring2020/Simulation_Partial/")

sourceCpp("./Cpp_R_Implementation/clickingAlgo.cpp")
source("./shared/allFunctions.R")
source("./shared/pseudoSamplingFuncs.R")

##########generate a dataset using Mallows distribution ###############
n<-20  
N<-20
rho0<-1:n
alpha0<-2

source("./shared/pseudoSamplingFuncs.R")
if(n>50){
  load("./shared/Cdfootrule.RData")
  fitvec = estimate_partition_function(alpha_vector = seq(0.01,10,0.2), n_items = n, metric = "footrule", nmc = 2000,degree=10)
}

ranking_data<-sample_mallows(rho0 = rho0,alpha0 = alpha0, n_samples = N)

################### fit pseudo mallows########################
n_samples<-200
rho_samples <- RhoSamples(n_samples = n_samples,data = ranking_data, alpha = alpha0, sigma = 0)


################visualisation###################
heatMatML<- heatMap(rho_samples,1:n)
par(mai=c(1,1,0.65,1))
image(heatMatML,col=tim.colors(64*10),zlim=c(0,1),axes=F,cex.lab=2, main = "Mallows")
par(mai=c(1,1,0.65,1))
image.plot(heatMatML, zlim=c(0,1),legend.only=T,horizontal = F)
