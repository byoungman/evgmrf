#include "evgmrf_types.h"

// CHOLMOD error handler: ignore "not positive definite" (detected via info()),
// pass any other CHOLMOD message on to R as a warning.
extern "C" void evgmrf_cholmod_handler(int status, const char* file,
                                       int line, const char* message) {
  if (status == CHOLMOD_NOT_POSDEF) return;
  Rf_warning("CHOLMOD status %d: '%s' at file '%s', line %d",
             status, message, file, line);
}

// ---- per-solver configuration, applied once at creation ----
inline void configure(SimLLT&) {}                 // Eigen's own code: nothing to do
inline void configure(SupLLT& c) {
  c.cholmod().error_handler = evgmrf_cholmod_handler;
}

// ---- log-determinant, per solver type ----
inline double chol_logdet(SimLLT& c) {
  // read diag(L) directly, without copying L
  return 2.0 * c.matrixL().nestedExpression().diagonal().array().log().sum();
}
inline double chol_logdet(SupLLT& c) {
  return c.logDeterminant();
}

// ---- generic implementations ----
template <typename Solver>
Solver& checked(Rcpp::XPtr<Solver>& chol) {
  if (chol.get() == nullptr) Rcpp::stop("Invalid factorisation object: re-run the analyse step.");
  return *chol;
}

template <typename Solver>
Solver& checked_factorised(Rcpp::XPtr<Solver>& chol) {
  Solver& c = checked(chol);
  if (c.info() != Eigen::Success)
    Rcpp::stop("No valid factorisation: the last factorize() failed or was never run.");
  return c;
}

template <typename Solver>
Rcpp::XPtr<Solver> analyze_impl(const Eigen::SparseMatrix<double>& A) {
  if (A.rows() != A.cols()) Rcpp::stop("Matrix A must be square.");
  Solver* chol = new Solver();
  configure(*chol);                  // install handler before any CHOLMOD call
  chol->analyzePattern(A);
  return Rcpp::XPtr<Solver>(chol, true);  // R's GC frees it
}

template <typename Solver>
double factorize_impl(Rcpp::XPtr<Solver> chol, const Eigen::SparseMatrix<double>& A) {
  Solver& c = checked(chol);
  if (A.rows() != A.cols()) Rcpp::stop("Matrix A must be square.");
  c.factorize(A);
  if (c.info() != Eigen::Success) return NA_REAL;   // not positive definite
  return chol_logdet(c);
}

// Dense right-hand side: numeric vector -> vector, numeric matrix -> matrix.
// Columns of a matrix are treated as separate right-hand sides.
template <typename Solver>
Rcpp::NumericVector solve_dense_impl(Rcpp::XPtr<Solver> chol, Rcpp::NumericVector b) {
  Solver& c = checked_factorised(chol);

  const bool is_mat = b.hasAttribute("dim");
  const int n = is_mat ? Rf_nrows(b) : b.size();
  const int k = is_mat ? Rf_ncols(b) : 1;
  if (c.rows() != n) Rcpp::stop("Dimension mismatch between factorisation and b.");

  Eigen::Map<const Eigen::MatrixXd> B(b.begin(), n, k);   // no copy
  Eigen::MatrixXd Z = c.solve(B);

  Rcpp::NumericVector out(Z.data(), Z.data() + Z.size());
  if (is_mat) out.attr("dim") = Rcpp::Dimension(n, k);
  return out;
}

// Sparse right-hand side: dgCMatrix in, dgCMatrix out.
template <typename Solver>
Eigen::SparseMatrix<double> solve_sparse_impl(Rcpp::XPtr<Solver> chol,
                                              const Eigen::Map<Eigen::SparseMatrix<double>> B) {
  Solver& c = checked_factorised(chol);
  if (c.rows() != B.rows()) Rcpp::stop("Dimension mismatch between factorisation and B.");
  Eigen::SparseMatrix<double> Z = c.solve(B);
  return Z;
}

// ---- quadratic forms: diag(B' A^{-1} B) = colSums((L^{-1} P B)^2) ----
// Only the forward solve with L is needed (half the work of a full solve), and
// only a length-k vector is returned to R, not an n x k matrix.

// Access CHOLMOD's factor, which Eigen keeps protected in CholmodBase.
struct SupFactorAccess : SupLLT {
  static cholmod_factor* get(SupLLT& c) {
    auto member = &SupFactorAccess::m_cholmodFactor;
    return c.*member;
  }
};

inline Eigen::VectorXd quadform(SimLLT& c, const Eigen::Map<Eigen::MatrixXd>& B) {
  Eigen::MatrixXd Y;
  if (c.permutationP().size() > 0) Y = c.permutationP() * B;
  else Y = B;
  c.matrixL().solveInPlace(Y);                       // Y = L^{-1} P B
  return Y.colwise().squaredNorm().transpose();
}

inline Eigen::VectorXd quadform(SupLLT& c, Eigen::Map<Eigen::MatrixXd>& B) {
  cholmod_factor* L = SupFactorAccess::get(c);
  cholmod_common& cm = c.cholmod();
  cholmod_dense b = Eigen::viewAsCholmod(B);         // no copy
  cholmod_dense* pb = M_cholmod_solve(CHOLMOD_P, L, &b, &cm);   // P B
  if (pb == nullptr) Rcpp::stop("CHOLMOD permutation failed.");
  cholmod_dense* y = M_cholmod_solve(CHOLMOD_L, L, pb, &cm);    // L^{-1} P B
  M_cholmod_free_dense(&pb, &cm);
  if (y == nullptr) Rcpp::stop("CHOLMOD triangular solve failed.");

  Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>> Y(
      static_cast<double*>(y->x), y->nrow, y->ncol, Eigen::OuterStride<>(y->d));
  Eigen::VectorXd out = Y.colwise().squaredNorm().transpose();
  M_cholmod_free_dense(&y, &cm);
  return out;
}

template <typename Solver>
Eigen::VectorXd quadform_impl(Rcpp::XPtr<Solver> chol, Rcpp::NumericMatrix b) {
  Solver& c = checked_factorised(chol);
  if (c.rows() != b.nrow()) Rcpp::stop("Dimension mismatch between factorisation and B.");
  Eigen::Map<Eigen::MatrixXd> B(b.begin(), b.nrow(), b.ncol());
  return quadform(c, B);
}

// ---- half solve: x = P' L^{-T} z, so that z ~ N(0, I) gives x ~ N(0, A^{-1}) ----
// Use this for simulation, not solve(): solve() gives A^{-1} z, with covariance A^{-2}.

inline Eigen::MatrixXd solve_Lt(SimLLT& c, const Eigen::Map<Eigen::MatrixXd>& Z) {
  Eigen::MatrixXd Y = Z;
  c.matrixU().solveInPlace(Y);                       // Y = L^{-T} Z  (matrixU() is L')
  if (c.permutationPinv().size() > 0) return c.permutationPinv() * Y;   // P^{-1} = P'
  return Y;
}

inline Eigen::MatrixXd solve_Lt(SupLLT& c, Eigen::Map<Eigen::MatrixXd>& Z) {
  cholmod_factor* L = SupFactorAccess::get(c);
  cholmod_common& cm = c.cholmod();
  cholmod_dense z = Eigen::viewAsCholmod(Z);         // no copy
  cholmod_dense* y = M_cholmod_solve(CHOLMOD_Lt, L, &z, &cm);   // L^{-T} Z
  if (y == nullptr) Rcpp::stop("CHOLMOD triangular solve failed.");
  cholmod_dense* x = M_cholmod_solve(CHOLMOD_Pt, L, y, &cm);    // P' L^{-T} Z
  M_cholmod_free_dense(&y, &cm);
  if (x == nullptr) Rcpp::stop("CHOLMOD permutation failed.");

  Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>> X(
      static_cast<double*>(x->x), x->nrow, x->ncol, Eigen::OuterStride<>(x->d));
  Eigen::MatrixXd out = X;
  M_cholmod_free_dense(&x, &cm);
  return out;
}

template <typename Solver>
Eigen::MatrixXd solve_Lt_impl(Rcpp::XPtr<Solver> chol, Rcpp::NumericMatrix z) {
  Solver& c = checked_factorised(chol);
  if (c.rows() != z.nrow()) Rcpp::stop("Dimension mismatch between factorisation and z.");
  Eigen::Map<Eigen::MatrixXd> Z(z.begin(), z.nrow(), z.ncol());
  return solve_Lt(c, Z);
}

// ---- exported: analyse ----
// [[Rcpp::export(.chol_analyze_simplicial)]]
Rcpp::XPtr<SimLLT> chol_analyze_simplicial(const Eigen::SparseMatrix<double>& A) {
  return analyze_impl<SimLLT>(A);
}

// [[Rcpp::export(.chol_analyze_supernodal)]]
Rcpp::XPtr<SupLLT> chol_analyze_supernodal(const Eigen::SparseMatrix<double>& A) {
  return analyze_impl<SupLLT>(A);
}

// ---- exported: simplicial ----
// [[Rcpp::export(.chol_factorize_simplicial)]]
double chol_factorize_simplicial(Rcpp::XPtr<SimLLT> chol, const Eigen::SparseMatrix<double>& A) {
  return factorize_impl(chol, A);
}

// [[Rcpp::export(.chol_solve_dense_simplicial)]]
Rcpp::NumericVector chol_solve_dense_simplicial(Rcpp::XPtr<SimLLT> chol, Rcpp::NumericVector b) {
  return solve_dense_impl(chol, b);
}

// [[Rcpp::export(.chol_solve_sparse_simplicial)]]
Eigen::SparseMatrix<double> chol_solve_sparse_simplicial(Rcpp::XPtr<SimLLT> chol,
                                                   const Eigen::Map<Eigen::SparseMatrix<double>> B) {
  return solve_sparse_impl(chol, B);
}

// [[Rcpp::export(.chol_quadform_simplicial)]]
Eigen::VectorXd chol_quadform_simplicial(Rcpp::XPtr<SimLLT> chol, Rcpp::NumericMatrix B) {
  return quadform_impl(chol, B);
}

// [[Rcpp::export(.chol_solve_Lt_simplicial)]]
Eigen::MatrixXd chol_solve_Lt_simplicial(Rcpp::XPtr<SimLLT> chol, Rcpp::NumericMatrix z) {
  return solve_Lt_impl(chol, z);
}

// [[Rcpp::export(.chol_L_simplicial)]]
Eigen::SparseMatrix<double> chol_L_simplicial(Rcpp::XPtr<SimLLT> chol) {
  SimLLT& c = checked_factorised(chol);
  return c.matrixL();
}

// ---- exported: supernodal ----
// [[Rcpp::export(.chol_factorize_supernodal)]]
double chol_factorize_supernodal(Rcpp::XPtr<SupLLT> chol, const Eigen::SparseMatrix<double>& A) {
  return factorize_impl(chol, A);
}

// [[Rcpp::export(.chol_solve_dense_supernodal)]]
Rcpp::NumericVector chol_solve_dense_supernodal(Rcpp::XPtr<SupLLT> chol, Rcpp::NumericVector b) {
  return solve_dense_impl(chol, b);
}

// [[Rcpp::export(.chol_solve_sparse_supernodal)]]
Eigen::SparseMatrix<double> chol_solve_sparse_supernodal(Rcpp::XPtr<SupLLT> chol,
                                                    const Eigen::Map<Eigen::SparseMatrix<double>> B) {
  return solve_sparse_impl(chol, B);
}

// [[Rcpp::export(.chol_quadform_supernodal)]]
Eigen::VectorXd chol_quadform_supernodal(Rcpp::XPtr<SupLLT> chol, Rcpp::NumericMatrix B) {
  return quadform_impl(chol, B);
}

// [[Rcpp::export(.chol_solve_Lt_supernodal)]]
Eigen::MatrixXd chol_solve_Lt_supernodal(Rcpp::XPtr<SupLLT> chol, Rcpp::NumericMatrix z) {
  return solve_Lt_impl(chol, z);
}

// ---- one-off log-determinant (no stored factorisation) ----
// [[Rcpp::export(.ldchol)]]
double ldchol(const Eigen::SparseMatrix<double>& A) {
  if (A.rows() != A.cols()) Rcpp::stop("Matrix A must be square.");
  SimLLT chol(A);
  if (chol.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed. Ensure A is positive definite.");
  return chol_logdet(chol);
}
