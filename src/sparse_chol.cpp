#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

// // Compute Cholesky decomposition and return the lower triangular factor L
// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> chol_sparse(const Eigen::SparseMatrix<double>& A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix A must be square.");
//   }
//   
//   // Compute Cholesky decomposition using Eigen
//   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
//   
//   if (chol.info() != Eigen::Success) {
//     Rcpp::stop("Cholesky decomposition failed. Ensure A is positive definite.");
//   }
//   
//   // Extract the lower triangular factor L
//   Eigen::SparseMatrix<double> L = chol.matrixL();
//   L.makeCompressed();  // Ensure it's in compressed format (CSC)
//   
//   return L;  // Return the lower triangular matrix L
// }

// Compute Cholesky decomposition, log(det(A)), and solve Az = b
// [[Rcpp::export(.cholAb)]]
Rcpp::List chol_logdet_solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) {
  if (A.rows() != A.cols()) {
    Rcpp::stop("Matrix A must be square.");
  }
  if (A.rows() != b.size()) {
    Rcpp::stop("Dimension mismatch: A and b must have compatible sizes.");
  }
  
  // Compute Cholesky decomposition using Eigen's SimplicialLLT
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
  
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed. Ensure A is positive definite.");
  }
  
  // Extract diagonal elements of L manually
  Eigen::SparseMatrix<double> L = chol.matrixL();
  Eigen::VectorXd diagL(L.rows());
  
  for (int i = 0; i < L.rows(); ++i) {
    diagL[i] = L.coeff(i, i);  // Extract diagonal value
  }
  
  // Compute log(det(A)) = 2 * sum(log(diagonal(L)))
  double logdet_A = 2.0 * diagL.array().log().sum();
  
  // Solve Az = b using the computed Cholesky decomposition
  Eigen::VectorXd z = chol.solve(b);
  
  // Return both logdet_A and solution z
  return Rcpp::List::create(
    Rcpp::Named("logdet_A") = logdet_A,
    Rcpp::Named("z") = z,
    Rcpp::Named("L") = L
  );
}

// Compute Cholesky decomposition, log(det(A)), and solve Az = b
// [[Rcpp::export(.ldchol)]]
double chol_logdet(const Eigen::SparseMatrix<double>& A) {
  if (A.rows() != A.cols()) {
    Rcpp::stop("Matrix A must be square.");
  }
  // Compute Cholesky decomposition using Eigen's SimplicialLLT
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
  
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed. Ensure A is positive definite.");
  }
  
  // Extract diagonal elements of L manually
  Eigen::SparseMatrix<double> L = chol.matrixL();
  Eigen::VectorXd diagL(L.rows());
  
  for (int i = 0; i < L.rows(); ++i) {
    diagL[i] = L.coeff(i, i);  // Extract diagonal value
  }
  
  // Compute log(det(A)) = 2 * sum(log(diagonal(L)))
  double logdet_A = 2.0 * diagL.array().log().sum();

  return logdet_A;
}

// // Solve a new system using an existing Cholesky factor L
// // [[Rcpp::export]]
// Eigen::VectorXd chol_solve_with_L(const Eigen::SparseMatrix<double>& L, const Eigen::VectorXd& b) {
//   if (L.rows() != L.cols()) {
//     Rcpp::stop("Matrix L must be square.");
//   }
//   if (L.rows() != b.size()) {
//     Rcpp::stop("Dimension mismatch: L and b must have compatible sizes.");
//   }
//   
//   // For a Cholesky factorization A = LL^T, to solve Ax = b:
//   // 1. Solve Ly = b for y using forward substitution
//   Eigen::VectorXd y = Eigen::VectorXd::Zero(b.size());
//   
//   for (int i = 0; i < L.rows(); i++) {
//     double sum = 0.0;
//     for (Eigen::SparseMatrix<double>::InnerIterator it(L, i); it; ++it) {
//       if (it.row() < i) {
//         sum += it.value() * y(it.row());
//       } else if (it.row() == i) {
//         y(i) = (b(i) - sum) / it.value();
//         break;
//       }
//     }
//   }
//   
//   // 2. Solve L^T x = y for x using backward substitution
//   Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
//   Eigen::SparseMatrix<double> Lt = L.transpose();
//   
//   for (int i = L.rows() - 1; i >= 0; i--) {
//     double sum = 0.0;
//     for (Eigen::SparseMatrix<double>::InnerIterator it(Lt, i); it; ++it) {
//       if (it.row() < i) {
//         sum += it.value() * x(it.row());
//       }
//     }
//     x(i) = (y(i) - sum) / L.coeff(i, i);
//   }
//   
//   return x;
// }
// 
// // Compute the Cholesky decomposition of a sparse matrix
// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> chol_sparse(Eigen::SparseMatrix<double> A) {
//   // Ensure the matrix is symmetric and positive definite
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
// 
//   // Perform Cholesky decomposition using Eigen's SimplicialLDLT
//     Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
// 
//     if (chol.info() != Eigen::Success) {
//         Rcpp::stop("Cholesky decomposition failed. Ensure the matrix is positive definite.");
//     }
// 
//     // Extract the lower triangular Cholesky factor
//     return chol.matrixL();
// }

// // Compute the Cholesky decomposition of a sparse matrix
// // [[Rcpp::export]]
// Eigen::VectorXd  chol_solve(Eigen::SparseMatrix<double> A, const Eigen::VectorXd& b) {
//   // Ensure the matrix is symmetric and positive definite
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   // Perform Cholesky decomposition using Eigen's SimplicialLDLT
//   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
//   
//   if (chol.info() != Eigen::Success) {
//     Rcpp::stop("Cholesky decomposition failed. Ensure the matrix is positive definite.");
//   }
//   
//   // Solve for x using the precomputed Cholesky decomposition
//   return chol.solve(b);
// }

// // Solve Ax = b using an existing Cholesky decomposition (SimplicialLLT)
// // [[Rcpp::export]]
// Eigen::VectorXd chol_solve_llt(const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol, Eigen::VectorXd b) {
//   // Solve for x using the existing Cholesky decomposition
//   return chol.solve(b);
// }

// // Define a class to store the Cholesky decomposition and allow repeated solving
// class CholeskySolver {
// private:
//   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol;  // Cholesky decomposition
//   
// public:
//   // Constructor: Compute Cholesky decomposition once
//   CholeskySolver(const Eigen::SparseMatrix<double>& A) {
//     if (A.rows() != A.cols()) {
//       Rcpp::stop("Matrix must be square.");
//     }
//     chol.compute(A);
//     if (chol.info() != Eigen::Success) {
//       Rcpp::stop("Cholesky decomposition failed. Ensure the matrix is positive definite.");
//     }
//   }
//   
//   // Solve Ax = b using the stored Cholesky decomposition
//   Eigen::VectorXd solve(const Eigen::VectorXd& b) {
//     return chol.solve(b);
//   }
// };
// 
// // Register the class as an Rcpp module
// RCPP_MODULE(CholeskySolverModule) {
//   Rcpp::class_<CholeskySolver>("CholeskySolver")
//   .constructor<Eigen::SparseMatrix<double>>()  // Compute Cholesky
//   .method("solve", &CholeskySolver::solve);    // Solve Ax = b
// }
// 
// // Compute Cholesky decomposition and return a dCHMsimpl object
// // [[Rcpp::export]]
// Rcpp::S4 chol_dCHMsimpl(Eigen::SparseMatrix<double> A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   // Perform Cholesky decomposition using Eigen
//   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
//   
//   if (chol.info() != Eigen::Success) {
//     Rcpp::stop("Cholesky decomposition failed. Ensure A is positive definite.");
//   }
//   
//   // Extract the L factor (lower triangular)
//   Eigen::SparseMatrix<double> L = chol.matrixL();
//   
//   // Convert Eigen sparse matrix to dgCMatrix for R compatibility
//   Rcpp::S4 dCHM("dCHMsimpl");
//   dCHM.slot("x") = Rcpp::NumericVector(L.valuePtr(), L.valuePtr() + L.nonZeros());
//   dCHM.slot("p") = Rcpp::IntegerVector(L.outerIndexPtr(), L.outerIndexPtr() + L.cols() + 1);
//   dCHM.slot("i") = Rcpp::IntegerVector(L.innerIndexPtr(), L.innerIndexPtr() + L.nonZeros());
//   dCHM.slot("perm") = Rcpp::IntegerVector::create();  // Empty permutation (no pivoting)
//   dCHM.slot("Dim") = Rcpp::IntegerVector::create(A.rows(), A.cols());
//   
//   return dCHM;
// }

// // Convert an Eigen SimplicialLLT Cholesky decomposition to dCHMsimpl for Matrix::solve()
// // [[Rcpp::export]]
// Rcpp::S4 convert_to_dCHMsimpl(const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol) {
//   // Extract the lower triangular factor L from existing Cholesky decomposition
//   Eigen::SparseMatrix<double> L = chol.matrixL();
//   L.makeCompressed();  // Ensure CSC format (required by Matrix package)
//   
//   // Create an S4 object of class dCHMsimpl
//   Rcpp::S4 dCHM("dCHMsimpl");
//   
//   // Assign correct slots for dCHMsimpl
//   dCHM.slot("x") = Rcpp::NumericVector(L.valuePtr(), L.valuePtr() + L.nonZeros());
//   dCHM.slot("p") = Rcpp::IntegerVector(L.outerIndexPtr(), L.outerIndexPtr() + L.cols() + 1);
//   dCHM.slot("i") = Rcpp::IntegerVector(L.innerIndexPtr(), L.innerIndexPtr() + L.nonZeros());
//   dCHM.slot("perm") = Rcpp::IntegerVector::create();  // No permutation
//   dCHM.slot("Dim") = Rcpp::IntegerVector::create(L.rows(), L.cols());
//   
//   return dCHM;
// }

// // Compute Cholesky decomposition and return as an RcppEigen object
// // [[Rcpp::export]]
// Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol_decompose(const Eigen::SparseMatrix<double>& A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   // Compute Cholesky decomposition using Eigen
//   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
//   
//   if (chol.info() != Eigen::Success) {
//     Rcpp::stop("Cholesky decomposition failed. Ensure A is positive definite.");
//   }
//   
//   return chol;
// }
// 
// // Solve Ax = b using the precomputed Cholesky decomposition
// // [[Rcpp::export]]
// Eigen::VectorXd chol_solve(const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol, const Eigen::VectorXd& b) {
//   if (chol.rows() != b.size()) {
//     Rcpp::stop("Dimension mismatch: A must have the same number of rows as b.");
//   }
//   
//   return chol.solve(b);
// }

// // Solve Ax = b using the precomputed Cholesky decomposition (passed as L)
// // [[Rcpp::export]]
// Eigen::VectorXd chol_solve(const Eigen::SparseMatrix<double>& L, const Eigen::VectorXd& b) {
//   if (L.rows() != b.size()) {
//     Rcpp::stop("Dimension mismatch: A must have the same number of rows as b.");
//   }
//   
//   // Solve for x using the lower triangular matrix L (Cholesky)
//   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(L);  // Recreate the solver with L
//   return chol.solve(b);
// }
// 
// // Compute Cholesky decomposition manually using a sparse matrix
// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> manual_chol_sparse(const Eigen::SparseMatrix<double>& A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   int n = A.rows();
//   Eigen::SparseMatrix<double> L(n, n);  // Initialize lower triangular matrix
//   std::vector<Eigen::Triplet<double>> tripletList;
//   
//   for (int j = 0; j < n; ++j) {
//     double sum_diag = 0.0;
//     
//     // Compute diagonal element
//     for (Eigen::SparseMatrix<double>::InnerIterator it(L, j); it; ++it) {
//       sum_diag += it.value() * it.value();
//     }
//     
//     double diag_val = std::sqrt(A.coeff(j, j) - sum_diag);
//     if (diag_val <= 0.0) {
//       Rcpp::stop("Matrix is not positive definite.");
//     }
//     
//     tripletList.push_back(Eigen::Triplet<double>(j, j, diag_val));
//     
//     // Compute lower triangular values
//     for (int i = j + 1; i < n; ++i) {
//       if (A.coeff(i, j) != 0) {  // Only compute nonzero elements
//         double sum_off_diag = 0.0;
//         
//         for (int k = 0; k < j; ++k) {
//           sum_off_diag += L.coeff(i, k) * L.coeff(j, k);
//         }
//         
//         double val = (A.coeff(i, j) - sum_off_diag) / diag_val;
//         tripletList.push_back(Eigen::Triplet<double>(i, j, val));
//       }
//     }
//   }
//   
//   // Convert triplets to a sparse matrix
//   L.setFromTriplets(tripletList.begin(), tripletList.end());
//   
//   return L;  // Return lower triangular matrix
// }

// // Manually compute Cholesky decomposition for a sparse matrix
// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> fixed_manual_chol_sparse(const Eigen::SparseMatrix<double>& A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   int n = A.rows();
//   std::vector<Eigen::Triplet<double>> tripletList;
//   
//   for (int j = 0; j < n; ++j) {
//     double sum_diag = 0.0;
//     
//     // Compute L(j, j) (diagonal elements)
//     for (int k = 0; k < j; ++k) {
//       double Ljk = 0.0;
//       for (auto it = tripletList.begin(); it != tripletList.end(); ++it) {
//         if (it->row() == j && it->col() == k) {
//           Ljk = it->value();
//           break;
//         }
//       }
//       sum_diag += Ljk * Ljk;
//     }
//     
//     double diag_val = std::sqrt(A.coeff(j, j) - sum_diag);
//     if (diag_val <= 0.0) {
//       Rcpp::stop("Matrix is not positive definite.");
//     }
//     tripletList.push_back(Eigen::Triplet<double>(j, j, diag_val));
//     
//     // Compute off-diagonal elements L(i, j)
//     for (int i = j + 1; i < n; ++i) {
//       if (A.coeff(i, j) != 0) {  // Only process nonzero elements
//         double sum_off_diag = 0.0;
//         
//         for (int k = 0; k < j; ++k) {
//           double Lik = 0.0, Ljk = 0.0;
//           
//           for (auto it = tripletList.begin(); it != tripletList.end(); ++it) {
//             if (it->row() == i && it->col() == k) {
//               Lik = it->value();
//             }
//             if (it->row() == j && it->col() == k) {
//               Ljk = it->value();
//             }
//           }
//           sum_off_diag += Lik * Ljk;
//         }
//         
//         double val = (A.coeff(i, j) - sum_off_diag) / diag_val;
//         tripletList.push_back(Eigen::Triplet<double>(i, j, val));
//       }
//     }
//   }
//   
//   // Convert triplet list into a sparse matrix
//   Eigen::SparseMatrix<double> L(n, n);
//   L.setFromTriplets(tripletList.begin(), tripletList.end());
//   
//   return L;
// }

// #include <omp.h>  // OpenMP for parallelization
// 
// // [[Rcpp::plugins(openmp)]]
// 
// // Compute Cholesky decomposition, log(det(A)), and solve Az = b in parallel
// // [[Rcpp::export]]
// Rcpp::List chol_logdet_solve_fully_parallel(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, int num_threads = 4) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix A must be square.");
//   }
//   if (A.rows() != b.size()) {
//     Rcpp::stop("Dimension mismatch: A and b must have compatible sizes.");
//   }
//   
//   // Set number of threads for Eigen
//   Eigen::setNbThreads(num_threads);
//   
//   // Compute Cholesky decomposition using Eigen's SimplicialLLT (Parallel in Eigen)
//   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol(A);
//   
//   if (chol.info() != Eigen::Success) {
//     Rcpp::stop("Cholesky decomposition failed. Ensure A is positive definite.");
//   }
//   
//   // Extract diagonal elements of L in parallel
//   Eigen::SparseMatrix<double> L = chol.matrixL();
//   Eigen::VectorXd diagL(L.rows());
//   
// #pragma omp parallel for num_threads(num_threads)
//   for (int i = 0; i < L.rows(); ++i) {
//     diagL[i] = L.coeff(i, i);  // Extract diagonal value
//   }
//   
//   // Compute log(det(A)) = 2 * sum(log(diagonal(L)))
//   double logdet_A = 2.0 * diagL.array().log().sum();
//   
//   // Solve Az = b in parallel using Cholesky decomposition
//   Eigen::VectorXd z;
// #pragma omp parallel
// {
// #pragma omp single nowait
//   z = chol.solve(b);
// }
// 
// // Return both logdet_A and solution z
// return Rcpp::List::create(
//   Rcpp::Named("logdet_A") = logdet_A,
//   Rcpp::Named("z") = z
// );
// }

// // Solve Ax = b using SuperLU
// // [[Rcpp::export]]
// Rcpp::List chol_logdet_solve_superlu(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, int num_threads = 4) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix A must be square.");
//   }
//   if (A.rows() != b.size()) {
//     Rcpp::stop("Dimension mismatch: A and b must have compatible sizes.");
//   }
//   
//   // Set number of threads for Eigen (SuperLU does not use OpenMP)
//   Eigen::setNbThreads(num_threads);
//   
//   // Compute LU decomposition using Eigen's SparseLU (SuperLU backend)
//   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//   solver.compute(A);
//   
//   if (solver.info() != Eigen::Success) {
//     Rcpp::stop("SuperLU decomposition failed.");
//   }
//   
//   // Compute log(det(A)) from LU decomposition (sum log U diagonal)
//   Eigen::VectorXd diagU = A.diagonal();
//   double logdet_A = diagU.array().log().sum();
//   
//   // Solve Az = b
//   Eigen::VectorXd z = solver.solve(b);
//   
//   return Rcpp::List::create(
//     Rcpp::Named("logdet_A") = logdet_A,
//     Rcpp::Named("z") = z
//   );
// }

// // Manually compute the SimplicialLLT sparse Cholesky decomposition
// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> manual_simplicialLLT(const Eigen::SparseMatrix<double>& A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   int n = A.rows();
//   Eigen::SparseMatrix<double> L(n, n);
//   std::vector<Eigen::Triplet<double>> tripletList;
//   
//   // Convert A to CSC format for efficient access
//   Eigen::SparseMatrix<double> A_csc = A;
//   A_csc.makeCompressed();
//   
//   // Storage for L (lower triangular)
//   std::vector<std::vector<std::pair<int, double>>> L_entries(n);
//   
//   for (int j = 0; j < n; ++j) {
//     double sum_diag = 0.0;
//     
//     // Compute L(j, j)
//     for (const auto& entry : L_entries[j]) {
//       sum_diag += entry.second * entry.second;
//     }
//     
//     double diag_val = std::sqrt(A_csc.coeff(j, j) - sum_diag);
//     if (diag_val <= 0.0) {
//       Rcpp::stop("Matrix is not positive definite.");
//     }
//     tripletList.push_back(Eigen::Triplet<double>(j, j, diag_val));
//     
//     // Compute lower triangular values
//     for (Eigen::SparseMatrix<double>::InnerIterator it(A_csc, j); it; ++it) {
//       int i = it.row();
//       if (i > j) {  // Only update lower triangular part
//         double sum_off_diag = 0.0;
//         
//         // Compute L(i, j)
//         for (const auto& entry : L_entries[j]) {
//           int k = entry.first;
//           double Ljk = entry.second;
//           double Lik = 0.0;
//           
//           // Find L(i, k) if it exists
//           for (const auto& e : L_entries[i]) {
//             if (e.first == k) {
//               Lik = e.second;
//               break;
//             }
//           }
//           sum_off_diag += Lik * Ljk;
//         }
//         
//         double val = (A_csc.coeff(i, j) - sum_off_diag) / diag_val;
//         tripletList.push_back(Eigen::Triplet<double>(i, j, val));
//         L_entries[i].push_back({j, val});
//       }
//     }
//   }
//   
//   // Convert triplets into a sparse matrix
//   L.setFromTriplets(tripletList.begin(), tripletList.end());
//   return L;
// }

// // Manually compute SimplicialLLT sparse Cholesky decomposition with proper error handling
// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> manual_simplicialLLT(const Eigen::SparseMatrix<double>& A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   int n = A.rows();
//   Eigen::SparseMatrix<double> L(n, n);
//   std::vector<Eigen::Triplet<double>> tripletList;
//   
//   // Convert A to CSC format for efficient access
//   Eigen::SparseMatrix<double> A_csc = A;
//   A_csc.makeCompressed();
//   
//   // Check if the matrix is symmetric (SimplicialLLT requires this)
//   if (!A_csc.isCompressed()) {
//     Rcpp::stop("Matrix is not in compressed form.");
//   }
//   
//   std::vector<std::vector<std::pair<int, double>>> L_entries(n);
//   
//   for (int j = 0; j < n; ++j) {
//     double sum_diag = 0.0;
//     
//     // Compute diagonal element L(j, j)
//     for (const auto& entry : L_entries[j]) {
//       sum_diag += entry.second * entry.second;
//     }
//     
//     double diag_val = A_csc.coeff(j, j) - sum_diag;
//     if (diag_val <= 0.0) {
//       Rcpp::stop("Matrix is not positive definite: zero or negative diagonal detected.");
//     }
//     diag_val = std::sqrt(diag_val);
//     tripletList.push_back(Eigen::Triplet<double>(j, j, diag_val));
//     
//     // Compute off-diagonal elements
//     for (Eigen::SparseMatrix<double>::InnerIterator it(A_csc, j); it; ++it) {
//       int i = it.row();
//       if (i > j) {  // Only update lower triangular part
//         double sum_off_diag = 0.0;
//         
//         // Compute L(i, j)
//         for (const auto& entry : L_entries[j]) {
//           int k = entry.first;
//           double Ljk = entry.second;
//           double Lik = 0.0;
//           
//           // Find L(i, k) if it exists
//           for (const auto& e : L_entries[i]) {
//             if (e.first == k) {
//               Lik = e.second;
//               break;
//             }
//           }
//           sum_off_diag += Lik * Ljk;
//         }
//         
//         double val = (A_csc.coeff(i, j) - sum_off_diag) / diag_val;
//         tripletList.push_back(Eigen::Triplet<double>(i, j, val));
//         L_entries[i].push_back({j, val});
//       }
//     }
//   }
//   
//   // Convert triplets into a sparse matrix
//   L.setFromTriplets(tripletList.begin(), tripletList.end());
//   return L;
// }

// // Manually compute SimplicialLLT sparse Cholesky decomposition with stopping criteria
// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> strict_simplicialLLT(const Eigen::SparseMatrix<double>& A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   int n = A.rows();
//   Eigen::SparseMatrix<double> L(n, n);
//   std::vector<Eigen::Triplet<double>> tripletList;
//   
//   // Convert A to CSC format for efficient access
//   Eigen::SparseMatrix<double> A_csc = A;
//   A_csc.makeCompressed();
//   
//   // Storage for L (lower triangular)
//   std::vector<std::vector<std::pair<int, double>>> L_entries(n);
//   
//   for (int j = 0; j < n; ++j) {
//     double sum_diag = 0.0;
//     
//     // Compute L(j, j)
//     for (const auto& entry : L_entries[j]) {
//       sum_diag += entry.second * entry.second;
//     }
//     
//     double diag_val = A_csc.coeff(j, j) - sum_diag;
//     
//     // 💡 If diag_val is zero or negative, STOP IMMEDIATELY (like Eigen::SimplicialLLT)
//     if (diag_val <= 0.0) {
//       Rcpp::stop("Matrix is not positive definite: found L(" + std::to_string(j) + ", " + std::to_string(j) + ") <= 0.");
//     }
//     
//     diag_val = std::sqrt(diag_val);
//     tripletList.push_back(Eigen::Triplet<double>(j, j, diag_val));
//     
//     // Compute off-diagonal elements
//     for (Eigen::SparseMatrix<double>::InnerIterator it(A_csc, j); it; ++it) {
//       int i = it.row();
//       if (i > j) {  // Only update lower triangular part
//         double sum_off_diag = 0.0;
//         
//         // Compute L(i, j)
//         for (const auto& entry : L_entries[j]) {
//           int k = entry.first;
//           double Ljk = entry.second;
//           double Lik = 0.0;
//           
//           // Find L(i, k) if it exists
//           for (const auto& e : L_entries[i]) {
//             if (e.first == k) {
//               Lik = e.second;
//               break;
//             }
//           }
//           sum_off_diag += Lik * Ljk;
//         }
//         
//         double val = (A_csc.coeff(i, j) - sum_off_diag) / diag_val;
//         tripletList.push_back(Eigen::Triplet<double>(i, j, val));
//         L_entries[i].push_back({j, val});
//       }
//     }
//   }
//   
//   // Convert triplets into a sparse matrix
//   L.setFromTriplets(tripletList.begin(), tripletList.end());
//   return L;
// }

// #include <cmath>
// 
// // Left-looking Sparse Cholesky Decomposition (Alternative to SimplicialLLT)
// // [[Rcpp::export]]
// Eigen::SparseMatrix<double> left_looking_cholesky(const Eigen::SparseMatrix<double>& A) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix must be square.");
//   }
//   
//   int n = A.rows();
//   Eigen::SparseMatrix<double> L(n, n);
//   std::vector<Eigen::Triplet<double>> tripletList;
//   
//   // Convert A to compressed column storage (CSC) for efficient access
//   Eigen::SparseMatrix<double> A_csc = A;
//   A_csc.makeCompressed();
//   
//   // Storage for L (lower triangular)
//   std::vector<std::vector<std::pair<int, double>>> L_entries(n);
//   
//   for (int j = 0; j < n; ++j) {
//     double sum_diag = 0.0;
//     
//     // Compute L(j, j)
//     for (const auto& entry : L_entries[j]) {
//       sum_diag += entry.second * entry.second;
//     }
//     
//     double diag_val = A_csc.coeff(j, j) - sum_diag;
//     
//     // If diag_val <= 0, the matrix is NOT positive definite → STOP
//     if (diag_val <= 0.0) {
//       Rcpp::stop("Matrix is not positive definite: found L(" + std::to_string(j) + ", " + std::to_string(j) + ") <= 0.");
//     }
//     
//     diag_val = std::sqrt(diag_val);
//     tripletList.push_back(Eigen::Triplet<double>(j, j, diag_val));
//     
//     // Compute off-diagonal elements (left-looking update)
//     for (Eigen::SparseMatrix<double>::InnerIterator it(A_csc, j); it; ++it) {
//       int i = it.row();
//       if (i > j) {  // Only update lower triangular part
//         double sum_off_diag = 0.0;
//         
//         // Compute L(i, j)
//         for (const auto& entry : L_entries[j]) {
//           int k = entry.first;
//           double Ljk = entry.second;
//           double Lik = 0.0;
//           
//           // Find L(i, k) if it exists
//           for (const auto& e : L_entries[i]) {
//             if (e.first == k) {
//               Lik = e.second;
//               break;
//             }
//           }
//           sum_off_diag += Lik * Ljk;
//         }
//         
//         double val = (A_csc.coeff(i, j) - sum_off_diag) / diag_val;
//         tripletList.push_back(Eigen::Triplet<double>(i, j, val));
//         L_entries[i].push_back({j, val});
//       }
//     }
//   }
//   
//   // Convert triplets into a sparse matrix
//   L.setFromTriplets(tripletList.begin(), tripletList.end());
//   return L;
// }
// 
// 
// class CholeskyCache {
// private:
//   Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> perm;
//   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol;
//   bool analyzed = false;
//   bool factorized = false;
//   int n;
//   
// public:
//   // Compute and store AMD permutation
//   void compute_permutation(const Eigen::SparseMatrix<double>& A) {
//     if (A.rows() != A.cols()) {
//       Rcpp::stop("Matrix must be square.");
//     }
//     n = A.rows();
//     Eigen::AMDOrdering<int> ordering;
//     ordering(A.selfadjointView<Eigen::Lower>(), perm);
//   }
//   
//   // Analyze pattern once (symbolic factorization)
//   void analyze(const Eigen::SparseMatrix<double>& A) {
//     if (perm.size() == 0) {
//       compute_permutation(A);
//     }
//     Eigen::SparseMatrix<double> A_perm = perm.transpose() * A * perm;
//     chol.analyzePattern(A_perm);
//     analyzed = true;
//   }
//   
//   // Factorize numerical values (reuse symbolic pattern)
//   void factorize(const Eigen::SparseMatrix<double>& A) {
//     if (!analyzed) analyze(A);
//     Eigen::SparseMatrix<double> A_perm = perm.transpose() * A * perm;
//     chol.factorize(A_perm);
//     if (chol.info() != Eigen::Success) {
//       Rcpp::stop("Cholesky factorization failed. Ensure matrix is positive definite.");
//     }
//     factorized = true;
//   }
//   
//   // Solve Ax = b using cached factorization
//   Eigen::VectorXd solve(const Eigen::VectorXd& b) {
//     if (!factorized) {
//       Rcpp::stop("Call factorize() before solve().");
//     }
//     Eigen::VectorXd b_perm = perm.transpose() * b;
//     Eigen::VectorXd x_perm = chol.solve(b_perm);
//     return perm * x_perm;
//   }
//   
//   // Compute log-determinant using cached factor
//   double logdet() {
//     if (!factorized) {
//       Rcpp::stop("Call factorize() before logdet().");
//     }
//     Eigen::SparseMatrix<double> L = chol.matrixL();
//     double logdet = 0.0;
//     for (int i = 0; i < L.rows(); ++i) {
//       logdet += std::log(L.coeff(i, i));
//     }
//     return 2.0 * logdet;
//   }
// };
// 
// RCPP_MODULE(CholeskyCacheModule) {
//   Rcpp::class_<CholeskyCache>("CholeskyCache")
//   .constructor()
//   .method("compute_permutation", &CholeskyCache::compute_permutation)
//   .method("analyze", &CholeskyCache::analyze)
//   .method("factorize", &CholeskyCache::factorize)
//   .method("solve", &CholeskyCache::solve)
//   .method("logdet", &CholeskyCache::logdet);
// }
// 
// // [[Rcpp::export]]
// Rcpp::List chol_logdet_solve_cached(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) {
//   if (A.rows() != A.cols()) {
//     Rcpp::stop("Matrix A must be square.");
//   }
//   if (A.rows() != b.size()) {
//     Rcpp::stop("Dimension mismatch: A and b must have compatible sizes.");
//   }
//   
//   // Create and use cache
//   CholeskyCache cache;
//   cache.analyze(A);      // symbolic (reuse if already called externally)
//   cache.factorize(A);    // numeric
//   
//   double logdet_A = cache.logdet();
//   Eigen::VectorXd z = cache.solve(b);
//   
//   return Rcpp::List::create(
//     Rcpp::Named("logdet_A") = logdet_A,
//     Rcpp::Named("z") = z
//   );
// }
