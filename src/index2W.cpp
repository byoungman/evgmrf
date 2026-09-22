// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Build a sparse adjacency matrix from grid indices
//'
//' Constructs a sparse binary adjacency matrix \code{W} where
//' \code{W[i, j] == 1} if points \code{i} and \code{j} are
//' rook-adjacent (differ by exactly 1 in one coordinate and 0 in the
//' other), and 0 otherwise.
//'
//' @param ind A numeric matrix with two columns giving the row/column
//'   (or x/y) coordinates of each point. Each row is one point.
//'
//' @return A sparse matrix (\code{Matrix::dgCMatrix}) of dimension
//'   \code{nrow(ind)} by \code{nrow(ind)}, symmetric, with 1s marking
//'   adjacent pairs and 0s elsewhere.
//'
//' @examples
//' ind <- as.matrix(expand.grid(1:5, 1:5))
//' W <- index2W(ind)
//'
//' @export
// [[Rcpp::export]]
arma::sp_mat index2W(const arma::mat& ind) {
 
int n = ind.n_rows;
arma::sp_mat W(n, n);   // sparse, initialized to all zeros

for (int i = 0; i < n; ++i) {
  for (int j = i + 1; j < n; ++j) {
    double d1 = std::abs(ind(i, 0) - ind(j, 0));
    double d2 = std::abs(ind(i, 1) - ind(j, 1));
     
    if ((d1 == 1.0 && d2 == 0.0) || (d1 == 0.0 && d2 == 1.0)) {
      W(i, j) = 1.0;
      W(j, i) = 1.0;
    }
  }
}

return W;

}