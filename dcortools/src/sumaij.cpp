#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
double SUMAIJBIJ(IntegerVector & IY, NumericVector & X, NumericVector & Y, NumericVector & XY, NumericVector & SX_X, NumericVector & SX_Y, NumericVector & SX_XY, NumericVector & SY_X, NumericVector & SY_Y, NumericVector & SY_XY)
{
  /* variable declaration */
  int k, pos, pos0;
  int n = Y.size();
  int L = ceil(std::log((double)(n)) / std::log(2.0));
  int P = floor(std::pow(2.0, (double)(L+1)));
  double GAMMA1, GAMMAX, GAMMAY, GAMMAXY;
  double RES = 0;
  NumericVector S(P), T(P), U(P), V(P);
  IntegerVector p(L+1);
  
  
  p[0] = 1;
  for (int ell = 1; ell <= L; ell++) {
    p[ell] = p[ell - 1] << 1;
  }
  GAMMA1 = n - 1 - 2 * IY[0] - 2; 
  GAMMAX = SX_X[n-1] - X[0] - 2* SY_X[0] - 2 * SX_X[0];
  GAMMAY = SX_Y[n-1] - Y[0] - 2* SY_Y[0] - 2 * SX_Y[0];
  GAMMAXY = SX_XY[n-1] - XY[0] - 2* SY_XY[0] - 2 * SX_XY[0];
  RES += X[0] * Y[0] * GAMMA1 + GAMMAXY - X[0] * GAMMAY - Y[0] * GAMMAX;
    
  for (int ii = 1; ii <= n - 1; ii++) {
    GAMMA1 = GAMMAX = GAMMAY = GAMMAXY=0;
    pos0 = 0;
    k = IY[ii - 1] - 1;
    for (int ell = 0; ell <= L - 1; ell++) {
      pos = pos0 + k;
      S[pos] +=  1;
      T[pos] +=  X[ii - 1];
      U[pos] +=  Y[ii - 1];
      V[pos] +=  XY[ii - 1];
      k = k  >> 1;
      pos0 += p[L-ell];
    }
    pos0 = 0;
    k = (IY[ii] - 1);
    for (int ell = 0; ell <= L - 1; ell++) {
        if (k % 2 == 1) {
        pos = pos0 + k - 1;
        GAMMA1 += S[pos]; 
        GAMMAX += T[pos];
        GAMMAY += U[pos];
        GAMMAXY+= V[pos];
      }
      k = (k >> 1);
      pos0 += p[L-ell];
    }
    GAMMA1 = 4 * GAMMA1 + n - 1 - 2 * IY[ii] - 2 * (ii + 1); 
    GAMMAX = 4 * GAMMAX + SX_X[n-1] - X[ii] - 2* SY_X[ii] - 2 * SX_X[ii];
    GAMMAY = 4 * GAMMAY + SX_Y[n-1] - Y[ii] - 2* SY_Y[ii] - 2 * SX_Y[ii];
    GAMMAXY= 4 * GAMMAXY + SX_XY[n-1] - XY[ii] - 2* SY_XY[ii] - 2 * SX_XY[ii];
    RES += X[ii] * Y[ii] * GAMMA1 + GAMMAXY - X[ii] * GAMMAY - Y[ii] * GAMMAX;
  }
  return RES;
}

