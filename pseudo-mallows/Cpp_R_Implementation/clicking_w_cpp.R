rm(list=ls())
require(RcppArmadillo)
require(Rcpp)
require(BayesMallows)
require(fields)
sourceCpp("./Cpp_R_Implementation/clickingAlgo.cpp")
source("./shared/allFunctions.R")
source("./shared/pseudoSamplingFuncs.R")
data_mode<-"simulation"

#################generate some data###################
if(data_mode=="simulation"){
  n<-50  
  N<-50
  rho0<-1:n
  alpha0<-3
  load("./shared/Cdfootrule.RData")
  if(n>20){
    fitvec = estimate_partition_function(alpha_vector = seq(0.01,10,0.2), n_items = 50,metric = "footrule", nmc = 2000,degree=10)
  }
  origin_data<-sample_mallows(rho0 = rho0,alpha0 = alpha0, n_samples = N)
  
}else if (data_mode == "potatoes"){
  origin_data<-potato_visual
  rho0<-potato_true_ranking
  N<-dim(origin_data)[1]
  n<-dim(origin_data)[2]
}

n<-dim(origin_data)[2]
N<-dim(origin_data)[1]
lambda<-5
clicking_data<-createClickData(origin_data,lambda)

centre_inferred<-n+1-rank(apply(clicking_data,2,sum),ties.method = "first")
data_init<-t(apply(clicking_data,1,generateInit_oneuser_2,centre = centre_inferred))

n_samples<-3000
start = proc.time()
results<-PseudoForClicks2(n_samples=n_samples,clicking_data=clicking_data, data_init = data_init, alpha =20,sigma=0)
durationCpp= proc.time()-start

userArray3D_cpp<-results$individuals[,,(n_samples/2+1):n_samples]
#userArray3D_cpp<-results$individuals[,,3000:6000]
K<-2
recsCpp<-vector()
truthCpp<-vector()

for(userj in 1:N){
  print(userj)
  heatMat_userjCpp<-heatMap(t(userArray3D_cpp[userj,,]),1:n)
  recsCpp<-rbind(recsCpp,recommendK(heatMat_userjCpp,K, clicking_data[userj,],type="click"))
  truthCpp<-rbind(truthCpp, realTopK(clicking_data[userj,], origin_data[userj,],K,type="click"))
}

correctOrWrongCpp<-recsCpp
for(userj in 1:N){
  correctOrWrongCpp[userj,]<-recsCpp[userj,]%in%truthCpp[userj,]
  
}

mean(correctOrWrongCpp)

heatMat_pseudo_Cpp<-heatMap(results$rhoMat,rho0)
par(mai=c(1,1,0.65,1))
image(heatMat_pseudo_Cpp,col=tim.colors(64*10),zlim=c(0,1),axes=F,cex.lab=2, main="pseudo-likelihood")
par(mai=c(1,1,0.65,1))
image.plot(heatMat_pseudo_Cpp, zlim=c(0,1),legend.only=T,horizontal = F)

