createPartialData<-function(data,remove_rate,method,K,lambda){
  if(missing(K)){
    K<-0
  }
  if(missing(lambda)){
    lambda<-0
  }
  if(missing(remove_rate)){
    remove_rate<-0
  }
  if(method=="random"){
    return(data*(runif(length(as.vector(data)))>remove_rate))  
  }else if(method == "topK"){
    return(data*(data<=K))
  }
 if(method == "poisson"){
   for(i in 1:dim(data)[1]){
     data[i,]<-data[i,]*(data[i,]<=rpois(1,lambda))
   }
   return(data)
 } 
}
ind_support<-function(vec){
  n<-length(vec)
  tmp<-1:n
  return(tmp[!1:n %in% vec])
}

containsZero<-function(vec){
  return(sum(vec == 0)>0)
}

zeroLocations<-function(vec){
  return(which(vec==0))
}


sampleRho<-function(data,sigma,alpha0){
  n<-dim(data)[2]
  tmpRank<-rank(apply(data,2,mean),ties.method="first")
  In<-rank(rnorm(n, mean=generateVOrderings(tmpRank),sd = sigma))
  support<-1:n
  rho<-rep(0,n)
  for(i in 1:n){
    i_curr<-which(In==i)
    dist<-sapply(support, oneDimfootrule, Rs=data[,i_curr])
    log_num<-(-alpha0/(n)*(dist)) - max(-alpha0/(n)*(dist)) 
    log_denom<- log(sum(exp(log_num)))
    probs<-exp((log_num-log_denom))
    rand<-runif(1)
    indOfCdf<-length(support)+1-sum(rand<=cdfForSamples(probs))
    rho[i_curr]<-support[indOfCdf]
    support<-setdiff(support,rho[i_curr])
  }
  
  return(rho)
}

sampleForOneUser<-function(userVec,rho_curr,alpha0){
  zeroLocats<-zeroLocations(userVec)
  supports<-ind_support(userVec)
  num_missing<-length(zeroLocats)
  user_j<-userVec
  for(j in sample(num_missing)){
    dist<-abs(supports-rho_curr[zeroLocats[j]])
    log_num<-(-alpha0/(n)*(dist)) - max(-alpha0/(n)*(dist)) 
    log_denom<- log(sum(exp(log_num)))
    probs<-exp((log_num-log_denom))
    rand<-runif(1)
    indOfCdf<-length(supports)+1-sum(rand<=cdfForSamples(probs))
    user_j[zeroLocats[j]]<-supports[indOfCdf]
    supports<-setdiff(supports,user_j[zeroLocats[j]]) 
  }
  return(user_j)
}

recommendK<-function(posteriorTable, K, partialVec, type){
  if(type=="partial"){
    ranked<-max(partialVec)  
  }else if(type == "click"){
    ranked<-sum(partialVec)
  }
  
  if(K == 1){
    posteriorProbs<-posteriorTable[,(ranked+1)]
  }
  if(K>1){
    posteriorProbs<-apply(posteriorTable[,((ranked+1):(ranked+K))],1,sum)  
  }
  return(which(posteriorProbs>=sort(posteriorProbs,decreasing = TRUE)[K])[1:K])
}

realTopK<-function(data_partial, data_orig,K,type){
  if(type=="partial"){
    ranked<-max(data_partial)  
  }else if(type=="click"){
    ranked<-sum(data_partial)
  }
  range<-(ranked+1):(ranked+K)
  return(which(data_orig%in%range))
}

createClickData<-function(data,lambda,method){
  if(missing(method)){
    method<-"poisson"
  }
  n<-dim(data)[2]
  N<-dim(data)[1]
  clicking<-data
  if(method == "poisson"){
    for(i in 1:N){
      clicking[i,]<-as.numeric(clicking[i,]<=min(max(1,rpois(1,lambda)),n))
    }
  }
  return(clicking)
}

sampleForOneUserClicks<-function(user_vec,rho_curr,alpha0){
  n<-length(user_vec)
  zeroLocats<-which(user_vec<1)
  nonzeroLocats<-which(user_vec >= 1)
  n1<-length(nonzeroLocats)
  n2<-length(zeroLocats)
  totalClicks<-length(nonzeroLocats)
  support1<-1:n1 #clicked
  support2<-1:n2 #unclicked
  rho_curr1<-rank(rho_curr[nonzeroLocats])
  rho_curr2<-rank(rho_curr[zeroLocats])
  user_j_clicked<-user_vec[nonzeroLocats]
  user_j_unclicked<-user_vec[zeroLocats]
  user_j<-user_vec
  for(j in sample(1:length(nonzeroLocats))){
    dist<-abs(support1-rho_curr1[j])
    log_num<-(-alpha0/(n1)*(dist)) - max(-alpha0/(n1)*(dist)) 
    log_denom<- log(sum(exp(log_num)))
    probs<-exp((log_num-log_denom))
    rand<-runif(1)
    indOfCdf<-length(support1)+1-sum(rand<=cdfForSamples(probs))
    user_j_clicked[j]<-support1[indOfCdf]
    support1<-setdiff(support1,user_j_clicked[j]) 
  }
  for(k in sample(1:length(zeroLocats))){
    dist<-abs(support2-rho_curr2[k])
    log_num<-(-alpha0/(n2)*(dist)) - max(-alpha0/(n2)*(dist)) 
    log_denom<- log(sum(exp(log_num)))
    probs<-exp((log_num-log_denom))
    rand<-runif(1)
    indOfCdf<-length(support2)+1-sum(rand<=cdfForSamples(probs))
    user_j_unclicked[k]<-support2[indOfCdf]
    support2<-setdiff(support2,user_j_unclicked[k]) 
  }
  user_j[nonzeroLocats]<-user_j_clicked
  user_j[zeroLocats]<-user_j_unclicked+ n1
  return(user_j)
}


sampleForOneUserClicks2<-function(user_vec,rho_curr,alpha0){
  n<-length(user_vec)
  zeroLocats<-which(user_vec<1)
  nonzeroLocats<-which(user_vec >= 1)
  totalClicks<-length(nonzeroLocats)
  support1<-1:totalClicks #clicked
  support2<-1:length(zeroLocats) #unclicked
  rho_curr1<-rank(rho_curr[nonzeroLocats])
  rho_curr2<-rank(rho_curr[zeroLocats])
  user_j_clicked<-user_vec[nonzeroLocats]
  user_j_unclicked<-user_vec[zeroLocats]
  user_j<-user_vec
  V<-generateVOrderings(user_vec)
  V1<-rank(V[nonzeroLocats])
  V2<-rank(V[zeroLocats])
  for(j in V1){
    dist<-abs(support1-rho_curr1[j])
    log_num<-(-alpha0/(n)*(dist)) - max(-alpha0/(n)*(dist)) 
    log_denom<- log(sum(exp(log_num)))
    probs<-exp((log_num-log_denom))
    rand<-runif(1)
    indOfCdf<-length(support1)+1-sum(rand<=cdfForSamples(probs))
    user_j_clicked[j]<-support1[indOfCdf]
    support1<-setdiff(support1,user_j_clicked[j]) 
  }
  for(k in V2){
    dist<-abs(support2-rho_curr2[k])
    log_num<-(-alpha0/(n)*(dist)) - max(-alpha0/(n)*(dist)) 
    log_denom<- log(sum(exp(log_num)))
    probs<-exp((log_num-log_denom))
    rand<-runif(1)
    indOfCdf<-length(support2)+1-sum(rand<=cdfForSamples(probs))
    user_j_unclicked[k]<-support2[indOfCdf]
    support2<-setdiff(support2,user_j_unclicked[k]) 
  }
  user_j[nonzeroLocats]<-user_j_clicked
  user_j[zeroLocats]<-user_j_unclicked+ totalClicks
  return(user_j)
}



generateInit_oneuser<-function(user_vec){
  n<-length(user_vec)
  zeroLocats<-which(user_vec<1)
  nonzeroLocats<-which(user_vec >= 1)
  user_vec[nonzeroLocats]<-sample(1:length(nonzeroLocats))
  user_vec[zeroLocats]<-sample(1:length(zeroLocats)+length(nonzeroLocats))
  return(user_vec)
}

generateInit_oneuser_2<-function(user_vec,centre){
  n<-length(user_vec)
  zeroLocats<-which(user_vec<1)
  nonzeroLocats<-which(user_vec >= 1)
  user_vec[nonzeroLocats]<-rank(centre[nonzeroLocats])
  user_vec[zeroLocats]<-rank(centre[zeroLocats])+length(nonzeroLocats)
  return(user_vec)
}
