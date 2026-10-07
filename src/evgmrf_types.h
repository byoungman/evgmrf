#ifndef EVGMRF_TYPES_H
#define EVGMRF_TYPES_H

#include <RcppArmadillo.h>      // must come before anything that includes Rcpp.h
#include <RcppEigen.h>
#include <RcppEigenCholmod.h>
#include <Eigen/CholmodSupport>

typedef Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> SupLLT;
typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> SimLLT;

#endif
