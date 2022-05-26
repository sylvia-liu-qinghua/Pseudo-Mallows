#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std; 

rowvec SampleForRho(mat data, rowvec ordering_vec, double alpha);
rowvec Cdf(rowvec probs);
rowvec leapAndShift(rowvec rho, const int& L);

//import functions from R
Function ranking("rank");
Function rnorm_R("rnorm");

// [[Rcpp::export]]
rowvec SampleForRho(mat data, rowvec ordering_vec,double alpha){
  int n = data.n_cols;
  rowvec colsums = sum(data,0);
  rowvec ranked1= as<rowvec>(wrap(ranking(colsums, true,"first"))); // without gaussian noise
  //rowvec V_order = GenerateVOrdering(ranked1);
  //rowvec with_gaussian_noise = as<rowvec>(wrap(rnorm_R(n,V_order,sigma))); //importing normal function from R is not absolutely necessary
  //rowvec In = as<rowvec>(wrap(ranking(with_gaussian_noise, true,"first")));
  rowvec In = ordering_vec;
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

