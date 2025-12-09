// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include "../inst/include/dcortools_types.h"
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd squaresym(Eigen::MatrixXd A) {
  int n = A.cols();
  Eigen::MatrixXd B = Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().
  rankUpdate(A.adjoint());
  return(B);
}


// [[Rcpp::export]]
NumericMatrix normalize_matrix(NumericMatrix & M, NumericVector & cmM, double & mM) {
  NumericVector m(M.nrow(), mM);
  unsigned int ncol = M.ncol();
  for (unsigned int n = 0; n < ncol; ++n) {
    M.row(n) = M.row(n) - cmM;
    M.column(n) = M.column(n) - cmM + m;
  }
  return M;
}

// [[Rcpp::export]]
NumericMatrix hadamard_product(NumericMatrix & X, NumericMatrix & Y){
  unsigned int ncol = X.ncol();
  unsigned int counter1 = 0;
  unsigned int counter2 = 0;
  for (unsigned int j=0; j < (ncol * ncol); j++) {
    X[counter1++] *= Y[counter2++];
  }
  return X;
}


// [[Rcpp::export]]
double matrix_prod_sum2(NumericMatrix & X, NumericMatrix & Y){
  unsigned int ncol = X.ncol();
  int res = 0;
  for (unsigned int j=0; j < (ncol * ncol); j++) {
    res += X[j] * Y[j];
  }
  return res;
}



// [[Rcpp::export]]
double sum_hadamard_power2(NumericMatrix & X){
  double rres = 0;
  unsigned int n = X.nrow();
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = i + 1; j < n; j++) {
      rres += X(j, i) * X(j, i);
    }
  }
  return 2 * rres;
}

// [[Rcpp::export]]
double sum_hadamard_power3(NumericMatrix & X){
  double rres = 0;
  unsigned int n = X.nrow();
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = i + 1; j < n; j++) {
      rres += X(j, i) * X(j, i) * X(j, i);
    }
  }
  return 2 * rres;
}




// [[Rcpp::export]]
NumericVector vector_product(NumericVector & X, NumericVector & Y){
  unsigned int n = X.size();
  int counter = 0;
  for (unsigned int j = 0; j < n; j++) {
    X[counter++] *= Y[j];
  }
  return X;
}

// [[Rcpp::export]]
double matrix_prod_sum(NumericMatrix & X, NumericMatrix & Y){
  double rres = 0;
  unsigned int n = X.nrow();
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = i + 1; j < n; j++) {
      rres += X(j, i) * Y(j, i);
    }
  }
  return 2 * rres;
}

// [[Rcpp::export]]
double vector_prod_sum(NumericVector & X, NumericVector & Y){
  double rres = 0;
  unsigned int n = X.size();
  for (unsigned int i = 0; i < n; i++) {
    rres += X[i] * Y[i];
  }
  return rres;
}

// [[Rcpp::export]]
double specific_vector_prod_sum(NumericVector & X,
                                NumericVector & Y,
                                NumericVector & gamma_1,
                                NumericVector & gamma_X,
                                NumericVector & gamma_Y,
                                NumericVector & gamma_XY){
  double rres = 0;
  unsigned int n = X.size();
  for (unsigned int i = 0; i < n; i++) {
    rres += X[i] * Y[i] * gamma_1[i] + gamma_XY[i] - X[i] * gamma_Y[i] - Y[i] * gamma_X[i];
  }
  return rres;
}

// [[Rcpp::export]]
double matrix_sum(NumericMatrix & X) {
  double res = 0;
  int n = X.nrow();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      res += X(i, j);
    }
  }
  return 2 * res;
}

// [[Rcpp::export]]
double matrix_prod_sum_sample(const NumericMatrix X, const NumericMatrix Y, const IntegerVector s){
  unsigned int n = X.nrow();
  double res = 0;
  for (unsigned int j = 0; j < n; j++) {
    for (unsigned int i = (j + 1); i < n; i++) {
      res += X(i, j) * Y(s[i] - 1, s[j] - 1);
    }
  }
  return 2 * res;
}

// [[Rcpp::export]]
double matrix_prod_sum_subset(const NumericMatrix X, const NumericMatrix Y, const IntegerVector s){
  unsigned int n = s.size();
  double res = 0;
  for (unsigned int j = 0; j < n; j++) {
    for (unsigned int i = (j + 1); i < n; i++) {
      res += X(s[i] - 1, s[j] - 1) * Y(s[i] - 1, s[j] - 1);
    }
  }
  return 2 * res;
}

// [[Rcpp::export]]
NumericVector colsums_subset(const NumericMatrix X, const IntegerVector s){
  unsigned int n = s.size();
  NumericVector cm(n);
  for (unsigned int j = 0; j < n; j++) {
    cm[j] = 0;
    for (unsigned int i = 0; i < n; i++) {
      cm[j] += X(s[i] - 1, s[j] - 1);
    }
  }
  return cm;
}



// [[Rcpp::export]]
double vector_prod_sum_sample(const NumericVector X, const NumericVector Y, const IntegerVector s){
  unsigned int n = X.size();
  double res = 0;
  for (unsigned int j = 0; j < n; j++) {
    res += X[j] * Y[s[j] - 1];
  }
  return res;
}
