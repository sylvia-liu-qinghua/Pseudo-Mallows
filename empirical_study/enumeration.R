rm(list = ls())
library(gtools)
require(fields)
require(BayesMallows)
require(tidyr)
require(dplyr)
require(combinat)
require(Rcpp)
require(RcppArmadillo)
source("./shared/allFunctions.R")
source("./shared/pseudoSamplingFuncs.R")
sourceCpp("./cpp_for_SGD.cpp")

swap<-function(vec, position, L=1){
  n<-length(vec)
  if(vec[position]-L<=0){
    index<-which(vec == vec[position]+L)
    tmp<-vec[position]
    vec[position]<-vec[index]
    vec[index]<-tmp
  }else if(vec[position]+L>n){
    index<-which(vec == vec[position]-L)
    tmp<-vec[position]
    vec[position]<-vec[index]
    vec[index]<-tmp
  }else{
    if(runif(1)>0.5){
      index<-which(vec == vec[position]-L)
      tmp<-vec[position]
      vec[position]<-vec[index]
      vec[index]<-tmp
    }else{
      index<-which(vec == vec[position]+L)
      tmp<-vec[position]
      vec[position]<-vec[index]
      vec[index]<-tmp
    }
  }
  return(vec)
  
}

l_and_s_all_neighbours<-function(perm, pos){
  n<-length(perm)
  result_mat<-matrix(data=NA, nrow = (n-1), ncol=n)
  a<-perm[pos]
  all_pos<-setdiff(1:n, pos)
  i<-1
  for(target_pos in all_pos){
    result<-perm
    b<-perm[target_pos]
    c<-sign(a-b)
    if(c<0){
      indices<-which((perm>a) & (perm<=b))
      result[indices]<-perm[indices]-1
    }else{
      indices<-which((perm>=b) &(perm<a))
      result[indices]<-perm[indices]+1
    }
    result[pos]<-b
    result_mat[i,]<-result
    i<-i+1
  }
  return(result_mat)
  
}

l_and_s<-function(perm, pos, target_pos){
  n<-length(perm)
  a<-perm[pos]
  result<-perm
  b<-perm[target_pos]
  c<-sign(a-b)
  if(c<0){
    indices<-which((perm>a) & (perm<=b))
    result[indices]<-perm[indices]-1
  }else{
    indices<-which((perm>=b) &(perm<a))
    result[indices]<-perm[indices]+1
  }
  result[pos]<-b
  return(result)
}

l_s_with_neighbours<-function(perm,pos){
  n<-length(perm)
  result<-perm
  a<-perm[pos]
  if(a == n){
    b<-n-1
    target_pos<-which(perm == (n-1))
    result[pos]<-n-1
    result[target_pos]<-n
    
  }else{
    target<-a+1
    target_pos<-which(perm == target)
    result[pos]<-target
    result[target_pos]<-a
  }
  
  return(result)
  
}

perm_equal<-function(perm1,perm2){
  n<-length(perm1)
  return(sum(perm1==perm2)==n)
}

n<-5
allPermutations<-permutations(n,n)
alpha0<-1.5
rho0<-1:n
N<-100
full_data<-sample_mallows(rho0 = rho0,alpha0 = alpha0,n_samples = N)
nmc<-1000000
#############Mallows estimate, only need to run once #############
mallows_model=compute_mallows(rankings=full_data, alpha_init = alpha0, alpha_jump = nmc-1, nmc = nmc, rho_thinning=100)
long = mallows_model$rho
wide<-spread(long,item,value)
rho_samples<-wide[wide['iteration']>nmc/2, 3:(n+2)]

heatMatML<- heatMap(rho_samples,1:n)
par(mai=c(1,1,0.65,1))
image(heatMatML,col=tim.colors(64*10),zlim=c(0,1),axes=F,cex.lab=2, main = "Mallows")
par(mai=c(1,1,0.65,1))
image.plot(heatMatML, zlim=c(0,1),legend.only=T,horizontal = F)

probs_Mallows<-vector()
for(j in 1:dim(allPermutations)[1]){
  tmp_perm<-allPermutations[j,]
  probs_Mallows<-c(probs_Mallows,sum(apply(rho_samples,1,perm_equal,perm2 = tmp_perm))/dim(rho_samples)[1])
}

#########Go through every single combination, Record the KLs ##############
n_samples<-200
KLs_n_dim<-vector()
KLs_marginal<-vector()
for(iter in 1:dim(allPermutations)[1]){
  print(paste(iter, "/", factorial(n), "iterations completed"))
  probs_pseudo<-vector()
  rho_samples_Pseudo<-vector()
  ordering_current<-allPermutations[iter,]
  for(i in 1:n_samples){
    tmp = SampleForRho(full_data,ordering_current,alpha0)  
    rho_samples_Pseudo<-rbind(rho_samples_Pseudo,tmp)
  }
  heatMatPseudo<- heatMap(rho_samples_Pseudo,1:n)
  # for(j in 1:dim(allPermutations)[1]){
  #   tmp_perm<-allPermutations[j,]
  #   probs_pseudo<-c(probs_pseudo,sum(apply(rho_samples_Pseudo,1,perm_equal,perm2 = tmp_perm))/n_samples)
  # }
  # KLs_n_dim<-c(KLs_n_dim,KL_margin(probs_pseudo,probs_Mallows, margin = 0.0000000001))
  KLs_marginal<-c(KLs_marginal,KL_margin(heatMatPseudo,heatMatML, 0.0000000001))
}


KLs_marginal[apply(allPermutations,1,isV,1:n)]
summary(KLs_marginal)

