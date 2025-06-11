// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(.eig1)]]
double eig1(const arma::sp_mat& A) {
  arma::vec eigval;
  // Compute smallest eigenvalue
  bool success = arma::eigs_sym(eigval, A, 1, "sa");  // "sa" = smallest algebraic
  
  if (!success) {
    Rcpp::stop("eigs_sym did not converge.");
  }
  
  return eigval[0];
}