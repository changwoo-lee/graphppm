#include <Rcpp.h>
using namespace Rcpp;

// calculating psm efficiently

// [[Rcpp::export]]
IntegerMatrix uppertri_pcm(IntegerMatrix zsamples, int n, int nsamples){
  IntegerMatrix uppertri_pcm(n,n);
  for(int isamples=0; isamples<nsamples; isamples++){
    for(int i=0; i<n; i++){
      for(int j=(i+1) ; j<n; j++){ // n must be greater than 2
        if(zsamples(isamples,i) == zsamples(isamples,j)) uppertri_pcm(i, j)++;                 
      }            
    }
  }
  return(uppertri_pcm);
}

// [[Rcpp::export]]
void update_uppertri_pcm(IntegerMatrix uppertri_pcm, IntegerVector z, int n){
  for(int i=0; i<n; i++){
    for(int j=(i+1) ; j<n; j++){ // n must be greater than 2
      if(z(i) == z(j)) uppertri_pcm(i, j)++;                 
    }            
  }
}

// update upper trianglular log-psm 
// adding weight * 1(z_i = z_j) of (i,j)th component of psm using logsumexp trick

// [[Rcpp::export]]
void update_uppertri_logpsm(NumericMatrix uppertri_logpsm, IntegerVector z, int n, double logweight){

  for(int i=0; i<n; i++){
    for(int j=i ; j<n; j++){ // n must be greater than 2
      if(z(i) == z(j)){
        if(uppertri_logpsm(i,j) == R_NegInf){
          uppertri_logpsm(i,j) = logweight;
        }else{
          uppertri_logpsm(i, j) = uppertri_logpsm(i, j) + log(1 + exp(logweight - uppertri_logpsm(i, j)));
        }
      }
    }            
  }
}
// 
// 
// // [[Rcpp::export]]
// NumericMatrix update_uppertri_logpsm_copy(NumericMatrix &uppertri_logpsm, IntegerVector z, int n, double logweight){
//   for(int i=0; i<n; i++){
//     for(int j=i ; j<n; j++){ // n must be greater than 2
//       if(z(i) == z(j)){
//         if(uppertri_logpsm(i,j) == R_NegInf){
//           uppertri_logpsm(i,j) = logweight;
//         }else{
//           uppertri_logpsm(i, j) = uppertri_logpsm(i, j) + log(1 + exp(logweight - uppertri_logpsm(i, j)));
//         }
//       }
//     }            
//   }
//   return(uppertri_logpsm);
// }




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
// 
// /*** R
// cls3 <- rbind(c(1,1,2,2),c(1,1,2,2),c(1,2,2,2))
// 
// input = mcclust::comp.psm(cls3)*3
// input[lower.tri(input, diag =T)] = 0
// mode(input) = "integer"
// inputcopy = input
// input
// update_uppertri_pcm(input, c(2,2,1,1), 4)
// input
// input/4
// 
// # same as..
// cls4 <- rbind(c(1,1,2,2),c(1,1,2,2),c(1,2,2,2),c(2,2,1,1))
// uppertri_pcm(cls4, 4, 4)/4
// mcclust::comp.psm(cls4)
// 
// 
// loginput = log(mcclust::comp.psm(cls3)*3/4)
// exp(loginput)
// update_uppertri_logpsm(loginput, c(2,2,1,1), 4, log(1/4))
// exp(loginput)
// 
// */
