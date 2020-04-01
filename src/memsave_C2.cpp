#include <RcppArmadillo.h>
using namespace Rcpp;

double euclidean1(double Xi,  double Xj, double prm) {
  (void) prm;
  return (Xi - Xj) * (Xi - Xj);
}

double minkowski1(double Xi,  double Xj, double prm) {
  return pow(fabs(Xi - Xj), prm);
}

double gaussian2(double aij, double prm) {
  (void) prm;
  return 1 - exp(- 0.5 * aij);
}

double disc2(double aij, double prm) {
  (void) prm;
  int x = aij!=0;
  return x;
}

double euclidean2(double aij, double prm) {
  (void) prm;
  return sqrt(aij);
}

double boundsq2(double aij, double prm) {
  (void) prm;
  return aij / (1 + aij);
}

double gaussian2par(double aij, double prm) {
  (void) prm;
  return 1 - exp(- 0.5 * aij / prm / prm);
}

double boundsq2par(double aij, double prm) {
  (void) prm;
  return aij / (prm * prm + aij);
}

double minkowski2par(double aij, double prm) {
  double x = pow(aij, 1 / prm);
  return x;
}

double alpha2par(double aij, double prm) {
  return pow(aij, 0.5 * prm);
}

double absol(double aij,  double prm) {
  (void) prm;
  return fabs(aij);
}


double disc2vec(double aij, double prm) {
  (void) prm;
  int x = aij!=0;
  return x;
}


double gaussian2vec(double aij, double prm) {
  (void) prm;
  return 1 - exp(- 0.5 * aij * aij);
}

double boundsq2vec(double aij, double prm) {
  (void) prm;
  return (aij * aij) / (1 + aij * aij);
}

double gaussian2vecpar(double aij, double prm) {
  return 1 - exp(- 0.5 * aij * aij / prm / prm);
}

double boundsq2vecpar(double aij, double prm) {
  return (aij * aij) / (prm * prm + aij * aij);;
}

double alpha2vecpar(double aij, double prm) {
  return pow(fabs(aij), prm);
}


// [[Rcpp::export]]
List dcovtermsmemvec(NumericVector & X, NumericVector & Y, std::string & metrX, std::string & metrY, double prmX = 0, double prmY = 0, bool calcdcor = false, bool calcperm = false){
  unsigned int n = X.size();
  double aij = 0; 
  double bij =0;
  double Sa =0;
  double Sb =0;
  double Sab=0;
  double Saa=0;
  double Sbb=0;
  double aijbij=0;
  double aijaij=0;
  double bijbij=0;
  List res;
  double adotdot=0;
  double bdotdot=0;
  double (*fX2)(double,double);
  double (*fY2)(double,double);
  NumericVector aidot(n);
  NumericVector bidot(n);

  if (prmX == 0) {
    if (metrX == "gaussian") {
      fX2 = &gaussian2vec;
    } else if (metrX == "boundsq") {
      fX2 = &boundsq2vec;
    } else if (metrX == "discrete") {
      fX2 = &disc2vec;
    } else {
      fX2 = &absol;
    }
  } else {
    if (metrX == "gaussian") {
      fX2 = &gaussian2vecpar;
    } else if (metrX == "boundsq") {
      fX2 = &boundsq2vecpar;
    } else if (metrX == "alpha") {
      fX2 = &alpha2vecpar;
    } else if (metrX == "discrete") {
      fX2 = &disc2vec;
    } else {
      fX2 = &absol;
    }
  }
  
  if (prmY == 0) {
    if (metrY == "gaussian") {
      fY2 = &gaussian2vec;
    } else if (metrY == "boundsq") {
      fY2 = &boundsq2vec;
    } else if (metrY == "discrete") {
      fY2 = &disc2vec;
    } else {
      fY2 = &absol;
    }
  } else {
    if (metrY == "gaussian") {
      fY2 = &gaussian2vecpar;
    } else if (metrY == "boundsq") {
      fY2 = &boundsq2vecpar;
    } else if (metrY == "alpha") {
      fY2 = &alpha2vecpar;
    } else if (metrY == "discrete") {
      fY2 = &disc2vec;
    } else {
      fY2 = &absol;
    }
  }
  

  for (unsigned int i = 0; i < n; i++) {
    Sa = 0;
    Sb = 0;
    for (unsigned int j = 0; j < n; j++) {
      aij = fX2(X[i] - X[j],prmX); 
      bij = fY2(Y[i] - Y[j],prmY);
      Sa += aij;
      Sb += bij;
      aijbij += aij * bij;
      if (calcdcor) {
        aijaij += aij * aij;
        bijbij += bij * bij;
      }
    }
    
    if (calcdcor) {
      Saa += Sa*Sa;
      Sbb += Sb*Sb;
      if (calcperm) {
        aidot[i] = Sa;
        bidot[i] = Sb;
      }
    }
    
 
    
    Sab += Sa*Sb;
    adotdot += Sa;
    bdotdot += Sb;
  }
  
  res["aijbij"] = aijbij;
  res["Sab"] = Sab;
  
  if (calcdcor) {
    res["Saa"] = Saa;
    res["Sbb"] = Sbb;
    res["adotdot"] = adotdot;
    res["bdotdot"] = bdotdot;
    res["aijaij"] = aijaij;
    res["bijbij"] = bijbij;
    if (calcperm) {
      res["aidot"] = aidot;
      res["bidot"] = bidot;
    }
  } else {
    res["Tab"] = adotdot * bdotdot;
  }
  
  return(res);
}


// [[Rcpp::export]]
List dcovtermsmem(NumericMatrix & X, NumericMatrix & Y,  std::string & metrX, std::string & metrY, double prmX = 0, double prmY = 0, bool calcdcor = false, bool calcperm = false){
  unsigned int n = X.nrow();
  unsigned int p = X.ncol();
  unsigned int q = Y.ncol();
  double l2norm(NumericVector,int);
  double aij = 0;
  double bij =0;
  double Sa =0;
  double Sb =0;
  double Saa=0;
  double Sbb=0;
  double Sab=0;
  double aijbij=0;
  double aijaij=0;
  double bijbij=0;
  List res;
  double adotdot=0;
  double bdotdot=0;
  double (*fX1)(double, double, double);
  double (*fY1)(double, double, double);
  double (*fX2)(double, double);
  double (*fY2)(double, double);
  NumericVector aidot(n);
  NumericVector bidot(n);
  

  if (prmX == 0) {
    if (metrX == "gaussian") {
      fX1 = &euclidean1;
      fX2 = &gaussian2;
    } else if (metrX == "boundsq") {
      fX1 = &euclidean1;
      fX2 = &boundsq2;
    } else if (metrX == "discrete") {
      fX1 = &euclidean1;
      fX2 = &disc2;
    } else {
      fX1 = &euclidean1;
      fX2 = &euclidean2;
    }
  } else {
    if (metrX == "gaussian") {
      fX1 = &euclidean1;
      fX2 = &gaussian2par;
    } else if (metrX == "boundsq") {
      fX1 = &euclidean1;
      fX2 = &boundsq2par;
    } else if (metrX == "alpha") {
      fX1 = &euclidean1;
      fX2 = &alpha2par;
    } else if (metrX == "minkowski") {
      fX1 = &minkowski1;
      fX2 = &minkowski2par;
    } else if (metrX == "discrete") {
      fX1 = &euclidean1;
      fX2 = &disc2;
    } else {
      fX1 = &euclidean1;
      fX2 = &euclidean2;
    }
  }
  
  if (prmY == 0) {
    if (metrY == "gaussian") {
      fY1 = &euclidean1;
      fY2 = &gaussian2;
    } else if (metrY == "boundsq") {
      fY1 = &euclidean1;
      fY2 = &boundsq2;
    } else if (metrY == "discrete") {
      fY1 = &euclidean1;
      fY2 = &disc2;
    } else {
      fY1 = &euclidean1;
      fY2 = &euclidean2;
    }
  } else {
    if (metrY == "gaussian") {
      fY1 = &euclidean1;
      fY2 = &gaussian2par;
    } else if (metrY == "boundsq") {
      fY1 = &euclidean1;
      fY2 = &boundsq2par;
    } else if (metrY == "alpha") {
      fY1 = &euclidean1;
      fY2 = &alpha2par;
    } else if (metrY == "minkowski") {
      fY1 = &minkowski1;
      fY2 = &minkowski2par;
    } else if (metrY == "discrete") {
      fY1 = &euclidean1;
      fY2 = &disc2;
    } else {
      fY1 = &euclidean1;
      fY2 = &euclidean2;
    }
  }

  for (unsigned int i = 0; i < n; i++) {
    Sa = 0;
    Sb = 0;
    for (unsigned int j = 0; j < n; j++) {
      aij = 0;
      bij = 0;
      for (unsigned int k = 0; k < p; ++k) {
        aij += fX1(X(i,k), X(j,k), prmX);
      }
      for (unsigned int k = 0; k < q; ++k) {
        bij += fY1(Y(i,k), Y(j,k), prmY);
      }
      aij = fX2(aij, prmX);
      bij = fY2(bij, prmY);
      Sa += aij;
      Sb += bij;
      aijbij += aij * bij;
      if (calcdcor) {
        aijaij += aij * aij;
        bijbij += bij * bij;
      }
    }
    
    if (calcdcor) {
      Saa += Sa*Sa;
      Sbb += Sb*Sb;
      if (calcperm) {
        aidot[i] = Sa;
        bidot[i] = Sb;
      }
    }
    
    Sab += Sa*Sb;
    adotdot += Sa;
    bdotdot += Sb;
  }
  
  res["aijbij"] = aijbij;
  res["Sab"] = Sab;
  
  if (calcdcor) {
    res["Saa"] = Saa;
    res["Sbb"] = Sbb;
    res["adotdot"] = adotdot;
    res["bdotdot"] = bdotdot;
    res["aijaij"] = aijaij;
    res["bijbij"] = bijbij;
    if (calcperm) {
      res["aidot"] = aidot;
      res["bidot"] = bidot;
    }
  } else {
    res["Tab"] = adotdot * bdotdot;
  }
  
  return(res);
}




// [[Rcpp::export]]
List dvartermsmemvec(NumericVector & X, std::string & metrX, double prmX = 0) {
  unsigned int n = X.size();
  double aij = 0; 
  double Sa =0;
  double aijaij=0;
  List res;
  NumericVector aidot(n);
  double adotdot=0;
  double (*fX2)(double,double);

  if (prmX == 0) {
    if (metrX == "gaussian") {
      fX2 = &gaussian2vec;
    } else if (metrX == "boundsq") {
      fX2 = &boundsq2vec;
    } else if (metrX == "discrete") {
      fX2 = &disc2vec;
    } else {
      fX2 = &absol;
    }
  } else {
    if (metrX == "gaussian") {
      fX2 = &gaussian2vecpar;
    } else if (metrX == "boundsq") {
      fX2 = &boundsq2vecpar;
    } else if (metrX == "alpha") {
      fX2 = &alpha2vecpar;
    } else if (metrX == "discrete") {
      fX2 = &disc2vec;
    } else {
      fX2 = &absol;
    }
  }
  
 
  for (unsigned int i = 0; i < n; i++) {
    Sa = 0;
    for (unsigned int j = 0; j < n; j++) {
      aij = fX2(X[i] - X[j],prmX); 
      Sa += aij;
      aijaij += aij * aij;
    }
    
    aidot[i] = Sa;
    adotdot += Sa;
  }
  
  res["aijaij"] = aijaij;
  res["aidot"] = aidot;
  res["adotdot"] = adotdot;
  
  return(res);
}




// [[Rcpp::export]]
List dvartermsmem(NumericMatrix & X,  std::string & metrX, double prmX = 0) {
  unsigned int n = X.nrow();
  unsigned int p = X.ncol();
  double l2norm(NumericVector,int);
  double aij = 0;
  double Sa =0;
  double aijaij=0;
  NumericVector aidot(n);
  List res;
  double adotdot=0;
  double (*fX1)(double, double, double);
  double (*fX2)(double, double);

  
  if (prmX == 0) {
    if (metrX == "gaussian") {
      fX1 = &euclidean1;
      fX2 = &gaussian2;
    } else if (metrX == "boundsq") {
      fX1 = &euclidean1;
      fX2 = &boundsq2;
    } else if (metrX == "discrete") {
      fX1 = &euclidean1;
      fX2 = &disc2;
    } else {
      fX1 = &euclidean1;
      fX2 = &euclidean2;
    }
  } else {
    if (metrX == "gaussian") {
      fX1 = &euclidean1;
      fX2 = &gaussian2par;
    } else if (metrX == "boundsq") {
      fX1 = &euclidean1;
      fX2 = &boundsq2par;
    } else if (metrX == "alpha") {
      fX1 = &euclidean1;
      fX2 = &alpha2par;
    } else if (metrX == "minkowski") {
      fX1 = &minkowski1;
      fX2 = &minkowski2par;
    } else if (metrX == "discrete") {
      fX1 = &euclidean1;
      fX2 = &disc2;
    } else {
      fX1 = &euclidean1;
      fX2 = &euclidean2;
    }
  }
  
  for (unsigned int i = 0; i < n; i++) {
    Sa = 0;
    for (unsigned int j = 0; j < n; j++) {
      aij = 0;
       for (unsigned int k = 0; k< p; ++k) {
        aij += fX1(X(i,k), X(j,k), prmX);
      }
      aij = fX2(aij, prmX);
      Sa += aij;
      aijaij += aij * aij;
    }
    aidot[i] = Sa;
    adotdot += Sa;
  }
  
  res["aijaij"] = aijaij;
  res["aidot"] = aidot;
  res["adotdot"] = adotdot;
  
  return(res);
}






// [[Rcpp::export]]
double aijbijmemvec(NumericVector & X, NumericVector & Y, std::string & metrX, std::string & metrY, double prmX = 0, double prmY = 0) {
  unsigned int n = X.size();
  double aij = 0; 
  double bij =0;
  double aijbij=0;
  double (*fX2)(double,double);
  double (*fY2)(double,double);
  NumericVector aidot(n);
  NumericVector bidot(n);
  
  if (prmX == 0) {
    if (metrX == "gaussian") {
      fX2 = &gaussian2vec;
    } else if (metrX == "boundsq") {
      fX2 = &boundsq2vec;
    } else if (metrX == "discrete") {
      fX2 = &disc2vec;
    } else {
      fX2 = &absol;
    }
  } else {
    if (metrX == "gaussian") {
      fX2 = &gaussian2vecpar;
    } else if (metrX == "boundsq") {
      fX2 = &boundsq2vecpar;
    } else if (metrX == "alpha") {
      fX2 = &alpha2vecpar;
    } else if (metrX == "discrete") {
      fX2 = &disc2vec;
    } else {
      fX2 = &absol;
    }
  }
  
  if (prmY == 0) {
    if (metrY == "gaussian") {
      fY2 = &gaussian2vec;
    } else if (metrY == "boundsq") {
      fY2 = &boundsq2vec;
    } else if (metrY == "discrete") {
      fY2 = &disc2vec;
    } else {
      fY2 = &absol;
    }
  } else {
    if (metrY == "gaussian") {
      fY2 = &gaussian2vecpar;
    } else if (metrY == "boundsq") {
      fY2 = &boundsq2vecpar;
    } else if (metrY == "alpha") {
      fY2 = &alpha2vecpar;
    } else if (metrY == "discrete") {
      fY2 = &disc2vec;
    } else {
      fY2 = &absol;
    }
  }
  
  
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = (i+1); j < n; j++) {
      aij = fX2(X[i] - X[j],prmX); 
      bij = fY2(Y[i] - Y[j],prmY);
      aijbij += aij * bij;
    }  
 }

  return(2 * aijbij);
}



// [[Rcpp::export]]
double aijbijmem(NumericMatrix & X, NumericMatrix & Y,  std::string & metrX, std::string & metrY, double prmX = 0, double prmY = 0){
  unsigned int n = X.nrow();
  unsigned int p = X.ncol();
  unsigned int q = Y.ncol();
  double l2norm(NumericVector,int);
  double aij = 0;
  double bij =0;
  double aijbij=0;
  double (*fX1)(double, double, double);
  double (*fY1)(double, double, double);
  double (*fX2)(double, double);
  double (*fY2)(double, double);

  
  if (prmX == 0) {
    if (metrX == "gaussian") {
      fX1 = &euclidean1;
      fX2 = &gaussian2;
    } else if (metrX == "boundsq") {
      fX1 = &euclidean1;
      fX2 = &boundsq2;
    } else if (metrX == "discrete") {
      fX1 = &euclidean1;
      fX2 = &disc2;
    } else {
      fX1 = &euclidean1;
      fX2 = &euclidean2;
    }
  } else {
    if (metrX == "gaussian") {
      fX1 = &euclidean1;
      fX2 = &gaussian2par;
    } else if (metrX == "boundsq") {
      fX1 = &euclidean1;
      fX2 = &boundsq2par;
    } else if (metrX == "alpha") {
      fX1 = &euclidean1;
      fX2 = &alpha2par;
    } else if (metrX == "minkowski") {
      fX1 = &minkowski1;
      fX2 = &minkowski2par;
    } else if (metrX == "discrete") {
      fX1 = &euclidean1;
      fX2 = &disc2;
    } else {
      fX1 = &euclidean1;
      fX2 = &euclidean2;
    }
  }
  
  if (prmY == 0) {
    if (metrY == "gaussian") {
      fY1 = &euclidean1;
      fY2 = &gaussian2;
    } else if (metrY == "boundsq") {
      fY1 = &euclidean1;
      fY2 = &boundsq2;
    } else if (metrY == "discrete") {
      fY1 = &euclidean1;
      fY2 = &disc2;
    } else {
      fY1 = &euclidean1;
      fY2 = &euclidean2;
    }
  } else {
    if (metrY == "gaussian") {
      fY1 = &euclidean1;
      fY2 = &gaussian2par;
    } else if (metrY == "boundsq") {
      fY1 = &euclidean1;
      fY2 = &boundsq2par;
    } else if (metrY == "alpha") {
      fY1 = &euclidean1;
      fY2 = &alpha2par;
    } else if (metrY == "minkowski") {
      fY1 = &minkowski1;
      fY2 = &minkowski2par;
    } else if (metrY == "discrete") {
      fY1 = &euclidean1;
      fY2 = &disc2;
    } else {
      fY1 = &euclidean1;
      fY2 = &euclidean2;
    }
  }
  
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = (i+1); j < n; j++) {
      aij = 0;
      bij = 0;
      for (unsigned int k = 0; k < p; ++k) {
        aij += fX1(X(i,k), X(j,k), prmX);
      }
      for (unsigned int k = 0; k < q; ++k) {
        bij += fY1(Y(i,k), Y(j,k), prmY);
      }
      aij = fX2(aij, prmX);
      bij = fY2(bij, prmY);
      aijbij += aij * bij;
    }
  }  
    

  return(2* aijbij);
}
