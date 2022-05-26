#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std; 

rowvec GenerateVOrdering(rowvec r);
rowvec SampleForOneUserClicks(rowvec clicking, rowvec rho_curr, double alpha, double sigma);
rowvec SampleForOneUserClicksV(rowvec clicking, rowvec rho_curr, double alpha, double sigma);
rowvec SampleRho(mat data, double alpha,double sigma);
double PseudoLogZn(double alpha, rowvec Rj, rowvec rho);
rowvec Cdf(rowvec probs);
mat RhoSamples(int n_samples,mat data,double alpha, double sigma);
List PseudoForClicks(int n_samples,mat clicking_data, mat data_init, double alphaInit, double sigma);
//import functions from R
Function ranking("rank");
Function rnorm_R("rnorm");

rowvec GenerateVOrdering(rowvec r){
  int n = r.n_elem;
  rowvec V_ordering = r;
  rowvec distToCentre = r - (r.max()+1)/2;
  if (n % 2 == 1){
    //ordering[which(distToCentre==0)]<-1
    int i = as_scalar(find(distToCentre==0));
    V_ordering.col(i) = 1;
  }
  vec positiveDist = distToCentre(find(distToCentre>0));
  for(int i = 0; i<positiveDist.n_elem; i++){
    float ru=float(rand())/RAND_MAX ; //a random uniform between 0 and 1
    if(ru<0.5){
      V_ordering(as_scalar(find(distToCentre == positiveDist(i)))) = positiveDist(i)*2+1;
      V_ordering(as_scalar(find(distToCentre == (-1)*positiveDist(i)))) = positiveDist(i)*2;
      
    }else {
      V_ordering(as_scalar(find(distToCentre == positiveDist(i)))) = positiveDist(i)*2;
      V_ordering(as_scalar(find(distToCentre == (-1)*positiveDist(i)))) = positiveDist(i)*2+1;
    }
  }
  
  return V_ordering;
}

rowvec SampleRho(mat data, double alpha,double sigma){
  int n = data.n_cols;
  rowvec colsums = sum(data,0);
  rowvec ranked1= as<rowvec>(wrap(ranking(colsums, true,"first"))); // without gaussian noise
  rowvec V_order = GenerateVOrdering(ranked1);
  rowvec with_gaussian_noise = as<rowvec>(wrap(rnorm_R(n,V_order,sigma))); //importing normal function from R is not absolutely necessary
  rowvec In = as<rowvec>(wrap(ranking(with_gaussian_noise, true,"first")));
  rowvec rho(n);
  Row<int> support(n);
  for(int i=0; i<n; i++){
    support(i) = i+1;
  }

  for(int i=0; i<n; i++){
    rowvec dists(support.n_elem);
    int i_curr = as_scalar(find(In == (i+1)));
    for(int j=0; j<support.n_elem; j++){
      dists(j) = sum(abs(data.col(i_curr)-support(j)));
    }
    rowvec log_num = (-alpha/n)*dists - max(-alpha/(n)*(dists)) ;
    double log_denom = log(sum(exp(log_num)));
    rowvec probs= exp((log_num-log_denom));
    
    double rand_u=double(rand())/RAND_MAX ; //a random uniform between 0 and 1

    int indOfCdf = support.n_elem - sum(rand_u <= Cdf(probs));
    rho(i_curr) = support(indOfCdf);
    support.shed_col(as_scalar(find(support == rho(i_curr))));
  }
  return rho;

}

rowvec Cdf(rowvec probs){
  int length = probs.n_elem;
  rowvec probs_normalized = probs / sum(probs); //in case probs do not add to 1
  rowvec cdf_result(length);
  double sum_curr = 0;
  for(int i=0; i<length;i++){
    cdf_result(i) = (sum_curr+=probs_normalized(i));
  }
  return(cdf_result);
}

// [[Rcpp::export]]
mat RhoSamples(int n_samples,mat data,double alpha, double sigma){
  int n =data.n_cols;
  mat rhoMat (n_samples, n, fill::zeros);
  for(int i = 0; i<n_samples; i++){
    if(i%500 == 0){
      cout<<i<<"iteration completed"<<endl;
    }
    rhoMat.row(i) = SampleRho(data, alpha,sigma);
  }
  return rhoMat;
}

rowvec SampleForOneUserClicks(rowvec clicking, rowvec rho_curr, double alpha, double sigma){
  //int n = clicking.n_elem;
  uvec zeroLocats = find(clicking <1);
  uvec nonzeroLocats = find(clicking>=1);
  int n2=zeroLocats.n_elem;
  int n1=nonzeroLocats.n_elem;
  rowvec support1(n1);
  for(int i=0; i<n1; i++){
    support1(i) = i+1;
  }
  rowvec support2(n2);
  for(int i=0; i<n2; i++){
    support2(i) = i+1;
  }
  rowvec tmp1 = rho_curr.cols(nonzeroLocats);
  rowvec rho_curr1 = as<rowvec>(wrap(ranking(tmp1,true,"first")));
  rowvec tmp2 = rho_curr.cols(zeroLocats);
  rowvec rho_curr2 = as<rowvec>(wrap(ranking(tmp2,true,"first")));
  rowvec user_j_clicked = clicking.cols(nonzeroLocats);
  rowvec user_j_unclicked = clicking.cols(zeroLocats);
  rowvec user_j_new = clicking;
  rowvec curr_ordering1 = shuffle(support1);  //simply using a uniform ordering for sampling
  rowvec curr_ordering2 = shuffle(support2);
  for(int j = 0; j<n1; j++){
    int curr_index = curr_ordering1[j]-1;
    rowvec dists1 = abs(support1 - rho_curr1(curr_index));
    rowvec log_num = (-alpha/n1*(dists1)) - max(-alpha/(n1)*(dists1));
    double log_denom = log(sum(exp(log_num)));
    rowvec probs = exp((log_num-log_denom));
    double rand_u=double(rand())/RAND_MAX ; 
    int indOfCdf = support1.n_elem - sum(rand_u <= Cdf(probs));
    user_j_new(nonzeroLocats(curr_index))=support1(indOfCdf);
    support1.shed_col(as_scalar(find(support1 == user_j_new(nonzeroLocats(curr_index)))));
  }
   for(int j = 0; j<n2; j++){
     int curr_index = curr_ordering2[j]-1;
     rowvec dists2 = abs(support2 - rho_curr2(curr_index));
     rowvec log_num = (-alpha/(n2)*(dists2)) - max(-alpha/(n2)*(dists2));
     double log_denom = log(sum(exp(log_num)));
     rowvec probs = exp((log_num-log_denom));
     double rand_u=double(rand())/RAND_MAX ; ; 
     int indOfCdf = support2.n_elem - sum(rand_u <= Cdf(probs));
     user_j_new(zeroLocats(curr_index))=support2(indOfCdf);
     support2.shed_col(as_scalar(find(support2 == user_j_new(zeroLocats(curr_index)))));
   }
  user_j_new(zeroLocats)+=n1;
  return user_j_new;
}


// [[Rcpp::export]]
double PseudoLogZn(double alpha, rowvec Rj, rowvec rho){
  int n = Rj.n_elem;
  rowvec support(n) ;
  support = rho;
  
  double logSum = 0;
  for(int i=0; i<n; i++){
    logSum+=log(sum(exp((abs(support - rho(i)))*(-alpha/n))));
    support.shed_col(as_scalar(find(support == Rj(i))));
  }
  return logSum;
}

// [[Rcpp::export]]
rowvec SampleForOneUserClicksV(rowvec clicking, rowvec rho_curr, double alpha, double sigma){
  //int n = clicking.n_elem;
  rowvec rho_V =GenerateVOrdering(rho_curr);
  uvec zeroLocats = find(clicking <1);
  uvec nonzeroLocats = find(clicking>=1);
  vec tmptmp1 = rho_V(nonzeroLocats);
  vec tmptmp2 = rho_V(zeroLocats);
  rowvec curr_ordering1 = as<rowvec>(wrap(ranking(tmptmp1,true,"first")));
  rowvec curr_ordering2 = as<rowvec>(wrap(ranking(tmptmp2,true,"first")));
  
  //cout<<curr_ordering1<<endl;
  //cout<<curr_ordering2<<endl;
  
  int n2=zeroLocats.n_elem;
  int n1=nonzeroLocats.n_elem;
  rowvec support1(n1);
  for(int i=0; i<n1; i++){
    support1(i) = i+1;
  }
  rowvec support2(n2);
  for(int i=0; i<n2; i++){
    support2(i) = i+1;
  }
  rowvec tmp1 = rho_curr.cols(nonzeroLocats);
  rowvec rho_curr1 = as<rowvec>(wrap(ranking(tmp1,true,"first")));
  rowvec tmp2 = rho_curr.cols(zeroLocats);
  rowvec rho_curr2 = as<rowvec>(wrap(ranking(tmp2,true,"first")));
  rowvec user_j_clicked = clicking.cols(nonzeroLocats);
  rowvec user_j_unclicked = clicking.cols(zeroLocats);
  rowvec user_j_new = clicking;
  //rowvec curr_ordering1 = shuffle(support1);  //simply using a uniform ordering for sampling
  //rowvec curr_ordering2 = shuffle(support2);
  for(int j = 0; j<n1; j++){
    //int curr_index = curr_ordering1[j]-1;
    int curr_index = as_scalar(find(curr_ordering1==(j+1)));
    //int curr_index = as_scalar(find(curr_ordering1==j) );
    rowvec dists1 = abs(support1 - rho_curr1(curr_index));
    rowvec log_num = (-alpha/n1*(dists1)) - max(-alpha/(n1)*(dists1));
    double log_denom = log(sum(exp(log_num)));
    rowvec probs = exp((log_num-log_denom));
    double rand_u=double(rand())/RAND_MAX ; 
    int indOfCdf = support1.n_elem - sum(rand_u <= Cdf(probs));
    user_j_new(nonzeroLocats(curr_index))=support1(indOfCdf);
    support1.shed_col(as_scalar(find(support1 == user_j_new(nonzeroLocats(curr_index)))));
  }
  for(int j = 0; j<n2; j++){
    //int curr_index = curr_ordering2[j]-1;
    int curr_index = as_scalar(find(curr_ordering2==(j+1)));
    rowvec dists2 = abs(support2 - rho_curr2(curr_index));
    rowvec log_num = (-alpha/(n2)*(dists2)) - max(-alpha/(n2)*(dists2));
    double log_denom = log(sum(exp(log_num)));
    rowvec probs = exp((log_num-log_denom));
    double rand_u=double(rand())/RAND_MAX ; ; 
    int indOfCdf = support2.n_elem - sum(rand_u <= Cdf(probs));
    user_j_new(zeroLocats(curr_index))=support2(indOfCdf);
    support2.shed_col(as_scalar(find(support2 == user_j_new(zeroLocats(curr_index)))));
  }
  user_j_new(zeroLocats)+=n1;
  return user_j_new;
}

// [[Rcpp::export]]
rowvec SampleForOneUserClicksV2(rowvec clicking, rowvec rho_curr, double alpha, double sigma){
  //int n = clicking.n_elem;
  rowvec rho_V =GenerateVOrdering(rho_curr);
  uvec zeroLocats = find(clicking <1);
  uvec nonzeroLocats = find(clicking>=1);
  vec tmptmp1 = rho_V(nonzeroLocats);
  vec tmptmp2 = rho_V(zeroLocats);
  rowvec curr_ordering1 = as<rowvec>(wrap(ranking(tmptmp1,true,"first")));
  curr_ordering1 = GenerateVOrdering(curr_ordering1);
  rowvec curr_ordering2 = as<rowvec>(wrap(ranking(tmptmp2,true,"first")));
  curr_ordering2 = GenerateVOrdering(curr_ordering2);
  //cout<<curr_ordering1<<endl;
  //cout<<curr_ordering2<<endl;
  
  int n2=zeroLocats.n_elem;
  int n1=nonzeroLocats.n_elem;
  rowvec support1(n1);
  for(int i=0; i<n1; i++){
    support1(i) = i+1;
  }
  rowvec support2(n2);
  for(int i=0; i<n2; i++){
    support2(i) = i+1;
  }
  rowvec tmp1 = rho_curr.cols(nonzeroLocats);
  rowvec rho_curr1 = as<rowvec>(wrap(ranking(tmp1,true,"first")));
  rowvec tmp2 = rho_curr.cols(zeroLocats);
  rowvec rho_curr2 = as<rowvec>(wrap(ranking(tmp2,true,"first")));
  rowvec user_j_clicked = clicking.cols(nonzeroLocats);
  rowvec user_j_unclicked = clicking.cols(zeroLocats);
  rowvec user_j_new = clicking;
  //rowvec curr_ordering1 = shuffle(support1);  //simply using a uniform ordering for sampling
  //rowvec curr_ordering2 = shuffle(support2);
  for(int j = 0; j<n1; j++){
    //int curr_index = curr_ordering1[j]-1;
    int curr_index = as_scalar(find(curr_ordering1==(j+1)));
    //int curr_index = as_scalar(find(curr_ordering1==j) );
    rowvec dists1 = abs(support1 - rho_curr1(curr_index));
    rowvec log_num = (-alpha/n1*(dists1)) - max(-alpha/(n1)*(dists1));
    double log_denom = log(sum(exp(log_num)));
    rowvec probs = exp((log_num-log_denom));
    double rand_u=double(rand())/RAND_MAX ; 
    int indOfCdf = support1.n_elem - sum(rand_u <= Cdf(probs));
    user_j_new(nonzeroLocats(curr_index))=support1(indOfCdf);
    support1.shed_col(as_scalar(find(support1 == user_j_new(nonzeroLocats(curr_index)))));
  }
  for(int j = 0; j<n2; j++){
    //int curr_index = curr_ordering2[j]-1;
    int curr_index = as_scalar(find(curr_ordering2==(j+1)));
    rowvec dists2 = abs(support2 - rho_curr2(curr_index));
    rowvec log_num = (-alpha/(n2)*(dists2)) - max(-alpha/(n2)*(dists2));
    double log_denom = log(sum(exp(log_num)));
    rowvec probs = exp((log_num-log_denom));
    double rand_u=double(rand())/RAND_MAX ; ; 
    int indOfCdf = support2.n_elem - sum(rand_u <= Cdf(probs));
    user_j_new(zeroLocats(curr_index))=support2(indOfCdf);
    support2.shed_col(as_scalar(find(support2 == user_j_new(zeroLocats(curr_index)))));
  }
  user_j_new(zeroLocats)+=n1;
  return user_j_new;
}


// [[Rcpp::export]]
List PseudoForClicks2(int n_samples,mat clicking_data, mat data_init, double alphaInit, double sigma){
  int N = clicking_data.n_rows;
  int n = clicking_data.n_cols;
  double alphaPrev = alphaInit;
  double alphaProp = alphaInit;
  mat rhoMat (n_samples, n, fill::zeros);
  vec alphaResults (n_samples);
  cube individuals(N, n, n_samples, fill::zeros); //n_rows, n_cols, n_slices, fill_type
  rowvec colsums = sum(data_init,0);
  rowvec rho_curr = as<rowvec>(wrap(ranking(colsums, true,"first"))); //initialized a rho based on current data
  mat data_curr = data_init;
  for(int i=0; i<n_samples; i++){
    if((i+1)%1000 == 0){
      cout<<"iteration"<< i+1 <<"completed"<<endl;
    }
    alphaResults(i) = alphaPrev;
    rho_curr= SampleRho(data_curr, alphaPrev,sigma);
    rhoMat.row(i) = rho_curr;
    for(int j =0; j<N; j++){
      data_curr.row(j) = SampleForOneUserClicksV(clicking_data.row(j), rho_curr, alphaPrev,sigma);
    }
    individuals.slice(i) = data_curr;
  }
  
  return List::create(Named("rhoMat")=rhoMat,
                      Named("individuals")=individuals,
                      Named("alphas")=alphaResults);
}

//uniform
// [[Rcpp::export]]
List PseudoForClicks(int n_samples,mat clicking_data, mat data_init, double alphaInit, double sigma){
  int N = clicking_data.n_rows;
  int n = clicking_data.n_cols;
  double alphaPrev = alphaInit;
  double alphaProp = alphaInit;
  mat rhoMat (n_samples, n, fill::zeros);
  vec alphaResults (n_samples);
  cube individuals(N, n, n_samples, fill::zeros); //n_rows, n_cols, n_slices, fill_type
  rowvec colsums = sum(data_init,0);
  rowvec rho_curr = as<rowvec>(wrap(ranking(colsums, true,"first"))); //initialized a rho based on current data
  mat data_curr = data_init;
  for(int i=0; i<n_samples; i++){
    if((i+1)%1000 == 0){
      cout<<"iteration"<< i+1 <<"completed"<<endl;
    }
    alphaResults(i) = alphaPrev;
    rho_curr= SampleRho(data_curr, alphaPrev,sigma);
    rhoMat.row(i) = rho_curr;
    for(int j =0; j<N; j++){
      data_curr.row(j) = SampleForOneUserClicks(clicking_data.row(j), rho_curr, alphaPrev,sigma);
    }
    individuals.slice(i) = data_curr;
  }
  
  return List::create(Named("rhoMat")=rhoMat,
                      Named("individuals")=individuals,
                      Named("alphas")=alphaResults);
}

//double V
// [[Rcpp::export]]
List PseudoForClicks3(int n_samples,mat clicking_data, mat data_init, double alphaInit, double sigma){
  int N = clicking_data.n_rows;
  int n = clicking_data.n_cols;
  double alphaPrev = alphaInit;
  double alphaProp = alphaInit;
  mat rhoMat (n_samples, n, fill::zeros);
  vec alphaResults (n_samples);
  cube individuals(N, n, n_samples, fill::zeros); //n_rows, n_cols, n_slices, fill_type
  rowvec colsums = sum(data_init,0);
  rowvec rho_curr = as<rowvec>(wrap(ranking(colsums, true,"first"))); //initialized a rho based on current data
  mat data_curr = data_init;
  for(int i=0; i<n_samples; i++){
    if((i+1)%1000 == 0){
      cout<<"iteration"<< i+1 <<"completed"<<endl;
    }
    alphaResults(i) = alphaPrev;
    rho_curr= SampleRho(data_curr, alphaPrev,sigma);
    rhoMat.row(i) = rho_curr;
    for(int j =0; j<N; j++){
      data_curr.row(j) = SampleForOneUserClicksV2(clicking_data.row(j), rho_curr, alphaPrev,sigma);
    }
    individuals.slice(i) = data_curr;
  }
  
  return List::create(Named("rhoMat")=rhoMat,
                      Named("individuals")=individuals,
                      Named("alphas")=alphaResults);
}
