rm(list = ls())

require(BayesMallows)
require(combinat)
require(Rcpp)
require(RcppArmadillo)
require(fields)
library(RcppHungarian)
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

#########generate data ##############
n<-5
rho0<-1:n
N<-100
alpha0<-1.5
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


##########initialization - leap and shift 5 times###############
tmp<-rank(apply(full_data,2,mean))
ordering_current<-generateVOrderings(tmp)

while(isV(ordering_current, tmp)){
  rand_1<-sample(1:(n-1),1)
  ordering_current<-l_s_with_neighbours(ordering_current, rand_1)
}

isV(ordering_current,tmp)

init<-ordering_current
KL_tracker<-vector()
ordering_tracker<-vector()
equal_tracker<-0
cost_matrices<-list()

for(iter in 1:30){
  print(iter)
  rho_current<-vector()
  for(i in 1:200){
    tmp = SampleForRho(full_data,ordering_current,alpha0)  
    rho_current<-rbind(rho_current,tmp)
  }
  margin_dist_ML<-heatMap(rho_samples,1:n)
  margin_dist_Pseudo<-heatMap(rho_current,1:n)
  KL_curr<-KL_margin(margin_dist_Pseudo,margin_dist_ML, 0.0000000001)
  KL_tracker<-c(KL_tracker,KL_curr)
  cost_matrix<-matrix(data=KL_curr, nrow = n, ncol = n)
  
  for(i in 1:n){
    for(j in setdiff(1:n,i)){
      perturbed<-l_and_s(ordering_current, pos = i,target_pos = j) 
      tmp_rhos<-vector()
      for(k in 1:100){
        tmp = SampleForRho(full_data,perturbed,alpha0)  
        tmp_rhos<-rbind(tmp_rhos,tmp)
      }
      margin_dist_Pseudo<-heatMap(tmp_rhos,1:n)
      KL<-KL_margin(margin_dist_Pseudo,margin_dist_ML, 0.0000000001)
      cost_matrix[i,j]<-KL #-KL_curr
    }
  }
  
  cost_matrices[[iter]]<-cost_matrix
  # cost2 <- cost_matrix-KL_curr
  soln <- HungarianSolver(cost_matrix)
  tmp<-soln$pairs
  # 
  # ##KL's difference
  # sum<-0
  # for(ind in 1:n){
  #   sum<-sum+cost2[tmp[ind,1],tmp[ind,2]]
  # }
  ordering_current<-ordering_current[soln$pairs[,2]]
  ordering_tracker<-rbind(ordering_tracker,ordering_current)
 
  if(iter>1){
    if(sum(ordering_tracker[iter-1,] == ordering_tracker[iter,])==n)  
      equal_tracker<-equal_tracker+1
    else{
      equal_tracker<-0
    }
    if(equal_tracker>=3){
      break
    }
  }
}

ordering_tracker<-rbind(init,ordering_tracker)

tmp<-rank(apply(full_data,2,mean))
ordering_tracker[,tmp]
