#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

/*
// [[Rcpp::export]]
NumericVector Simulation(const NumericMatrix grid,
                         const NumericMatrix Srow,
                         const NumericMatrix Scol,
                         const int nsim) {
  
  int nbeta = Srow.ncol();
  
  NumericVector ll(nsim);
  NumericVector llavg(grid.nrow());
  arma::vec mean(nbeta, fill::zeros);
  arma::mat Sl(nbeta,nbeta);
  arma::mat Slinv(nbeta,nbeta);

  for (int i=0; i < grid.nrow(); i++) {
    Sl = grid(i,1)*Srow + grid(i,2)*Scol;
    Slinv = arma::pinv(Sl);
    
    for (int j=0; j < nsim; j++) {
      arma::vec mvnrnd(mean, Slinv);
    }
  }
  
  return llavg;
}
 
*/

// [[Rcpp::export]]
NumericMatrix IndGreater(NumericVector x) {
  int n = x.size();
  NumericMatrix elem(n);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (x[j] >= x[i]) {
        elem(j,i) = 1;
      } else {
        elem(j,i) = 0;
      }
    }
  }
  return elem;
}

// [[Rcpp::export]]
NumericMatrix IndLess(NumericVector x) {
  int n = x.size();
  NumericMatrix elem(n);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (x[j] <= x[i]) {
        elem(j,i) = 1;
      } else {
        elem(j,i) = 0;
      }
    }
  }
  return elem;
}

// [[Rcpp::export]]
NumericMatrix IndEqual(NumericVector x) {
  int n = x.size();
  NumericMatrix elem(n);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (x[j] == x[i]) {
        elem(j,i) = 1;
      } else {
        elem(j,i) = 0;
      }
    }
  }
  return elem;
}

// [[Rcpp::export]]
int Ind2(NumericVector x, NumericVector y, double a, double b) {
  int n = x.size();
  int sum = 0;
  for (int i=0; i<n; i++) {
    if (x[i] >= a and y[i] >= b) {
      sum = sum + 1;
    } else {
      sum = sum + 0;
    }
  }
  return sum;
}

// [[Rcpp::export]]
NumericMatrix risksetC(NumericVector x, NumericVector y) {
  int n = x.size();
  NumericMatrix riskset(n);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      riskset (j,i) = Ind2(x, y, x[j], y[i]);
    }
  }
  return riskset;
}

// [[Rcpp::export]]
NumericMatrix DeltaC(NumericVector x, NumericVector y) {
  int n = x.size();
  NumericMatrix delta(n);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      delta (j,i) = x[j]*y[i];
    }
  }
  return delta;
}



// [[Rcpp::export]]
double logLikC(const NumericMatrix riskset,
               const NumericMatrix logtheta,
               const NumericMatrix delta,
               const NumericMatrix I1,
               const NumericMatrix I2,
               const NumericMatrix I3) {

  int n = riskset.nrow();
  double sum = 0;

  /* Calculation of L1 */
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (riskset(i,j) > 0) {
        sum = sum +
          delta(i,j)*I1(i,j)*(
              logtheta(i,j)*I3(i,j) -
              std::log(riskset(i,j) + I2(i,j)*(std::exp(logtheta(i,j)) - 1))
          );
      } else {sum = sum + 0;}
    }
  }
  
  /* To calculate L2:
   * I1 = t(I2)
   * I2 = t(I1)
   * I3 = I6 */
  return(-sum);
}


// // [[Rcpp::export]]
// double logLikC(const NumericMatrix riskset,
//                const NumericMatrix logtheta,
//                const NumericMatrix delta,
//                const NumericMatrix I1,
//                const NumericMatrix I2,
//                const NumericMatrix I3,
//                const NumericMatrix I4) {
//   
//   int n = riskset.nrow();
//   double sum1 = 0;
//   double sum2 = 0;
//   double result;
//   
//   /* Calculation of L1 */
//   for (int i=0; i<n; i++) {
//     for (int j=0; j<n; j++) {
//       if (riskset(j,i) > 0) {
//         sum1 = sum1 +
//           delta(j,i)*I1(i,j)*(
//               logtheta(j,i)*I3(i,j) -
//               std::log(riskset(j,i) + I2(i,j)*(std::exp(logtheta(j,i)) - 1))
//           );
//       } else {sum1 = sum1 + 0;}
//     } 
//   }
//   
//   /* Calculation of L2 */
//   /* Note that for L2 compared to L1 we have: I1 = t(I2), I2 = t(I1) */
//   for (int i=0; i<n; i++) {
//     for (int j=0; j<n; j++) {
//       if (riskset(i,j) > 0) {
//         sum2 = sum2 +
//           delta(i,j)*I2(j,i)*(
//               logtheta(i,j)*I4(i,j) -
//               std::log(riskset(i,j) + I1(j,i)*(std::exp(logtheta(i,j)) - 1))
//           );
//       } else {sum2 = sum2 + 0;}
//     } 
//   }
//   
//   result = -sum1 - sum2;
//   
//   return(result);
// }

// // [[Rcpp::export]]
// NumericVector gradientC(const NumericMatrix riskset,
//                         const NumericMatrix logtheta,
//                         const Rcpp::List deriv,
//                         const int df,
//                         const NumericMatrix delta,
//                         const NumericMatrix I1,
//                         const NumericMatrix I2,
//                         const NumericMatrix I3,
//                         const NumericMatrix I4) {
// 
//   int n = riskset.nrow();
//   int totalparam = df*df;
//   NumericVector result(totalparam);
// 
//   /* Transform list of derivative matrices into vector of matrices */
//   std::vector<NumericMatrix> deriv_vec(totalparam);
//   for(int k = 0; k < totalparam; ++k) {
//     NumericMatrix deriv_R = deriv[k];
//     /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
//      deriv_vec[k] = derivMat; */
//     deriv_vec[k] = deriv_R;
//   }
// 
//   for (int m=0; m<totalparam; m++) {
// 
//     double sum1 = 0;
// 
//     /* Calculation of L1 */
//     for (int i=0; i<n; i++) {
//       for (int j=0; j<n; j++) {
//         if (riskset(j,i) > 0) {
//           sum1 = sum1 +
//             delta(j,i)*I1(i,j)*deriv_vec[m](j,i)*(
//                 I3(i,j) -
//                 I2(i,j)*(std::exp(logtheta(j,i)))/(riskset(j,i) + I2(i,j)*(std::exp(logtheta(j,i)) - 1))
//             );
//         } else {sum1 = sum1 + 0;}
//       }
//     }
// 
//     result[m] = -sum1; /* Return derivative of -loglik */
// 
//   }
// 
//   return(result);
// 
// }
// 
// // [[Rcpp::export]]
// NumericMatrix hessianC(const NumericMatrix riskset,
//                        const NumericMatrix logtheta,
//                        const Rcpp::List deriv,
//                        const int df,
//                        const NumericMatrix delta,
//                        const NumericMatrix I1,
//                        const NumericMatrix I2) {
// 
//   int n = riskset.nrow();
//   int totalparam = df*df;
// 
//   NumericMatrix result(totalparam);
// 
//   /* Transform list of derivative matrices into vector of matrices */
//   std::vector<NumericMatrix> deriv_vec(totalparam);
//   for(int k = 0; k < totalparam; ++k) {
//     NumericMatrix deriv_R = deriv[k];
//     /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
//      deriv_vec[k] = derivMat; */
//     deriv_vec[k] = deriv_R;
//   }
// 
//   /* Calculate upper triangle of hessian */
//   for (int m = 0; m < totalparam; m++) {
//     for (int l = m; l < totalparam; l++) {
// 
//       double sum1 = 0;
// 
//       for (int i=0; i<n; i++) {
//         for (int j=0; j<n; j++) {
//           if (riskset(j,i) > 0) {
//             sum1 = sum1 +
//               delta(j,i)*I1(i,j)*deriv_vec[m](j,i)*deriv_vec[l](j,i)*(riskset(j,i) - I2(i,j))*I2(i,j)*std::exp(logtheta(j,i))/pow(riskset(j,i) - I2(i,j) + I2(i,j)*(std::exp(logtheta(j,i))),2);
//           } else {sum1 = sum1 + 0;}
//         }
//       }
// 
//       result(m,l) = sum1;
//       result(l,m) = result(m,l); /* Symmetric hessian: upper triangle = lower triangle */
// 
//     }
//   }
// 
//     return(result); /* Return second derivative of -loglik */
// }
  

// [[Rcpp::export]]
NumericVector gradientC(const NumericMatrix riskset,
                 const NumericMatrix logtheta,
                 const Rcpp::List deriv,
                 const int df,
                 const NumericMatrix delta,
                 const NumericMatrix I1,
                 const NumericMatrix I2,
                 const NumericMatrix I3,
                 const NumericMatrix I4) {

  int n = riskset.nrow();
  int totalparam = df*df;
  NumericVector result(totalparam);

  /* Transform list of derivative matrices into vector of matrices */
  std::vector<NumericMatrix> deriv_vec(totalparam);
  for(int k = 0; k < totalparam; ++k) {
    NumericMatrix deriv_R = deriv[k];
    /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
     deriv_vec[k] = derivMat; */
    deriv_vec[k] = deriv_R;
  }

  for (int m=0; m<totalparam; m++) {

    double sum1 = 0;
    double sum2 = 0;

    /* Calculation of L1 */
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (riskset(j,i) > 0) {
          sum1 = sum1 +
            delta(j,i)*I1(i,j)*deriv_vec[m](j,i)*(
                I3(i,j) -
                I2(i,j)*(std::exp(logtheta(j,i)))/(riskset(j,i) + I2(i,j)*(std::exp(logtheta(j,i)) - 1))
            );
        } else {sum1 = sum1 + 0;}
      }
    }

    /* Calculation of L2 */
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (riskset(i,j) > 0) {
          sum2 = sum2 +
            delta(i,j)*I2(j,i)*deriv_vec[m](i,j)*(
                I4(i,j) -
                I1(j,i)*(std::exp(logtheta(i,j)))/(riskset(i,j) + I1(j,i)*(std::exp(logtheta(i,j)) - 1))
            );
        } else {sum2 = sum2 + 0;}
      }
    }

    result[m] = -sum1-sum2;

  }

return(result);

}


// [[Rcpp::export]]
NumericVector gradientPoly(const NumericMatrix riskset,
                        const NumericMatrix logtheta,
                        const Rcpp::List deriv,
                        const int df,
                        const NumericMatrix delta,
                        const NumericMatrix I1,
                        const NumericMatrix I2,
                        const NumericMatrix I3,
                        const NumericMatrix I4) {
  
  int n = riskset.nrow();
  int totalparam = df;
  NumericVector result(totalparam);

  /* Transform list of derivative matrices into vector of matrices */
  std::vector<NumericMatrix> deriv_vec(totalparam);
  for(int k = 0; k < totalparam; ++k) {
    NumericMatrix deriv_R = deriv[k];
    /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
     deriv_vec[k] = derivMat; */
    deriv_vec[k] = deriv_R;
  }

  for (int m=0; m<totalparam; m++) {
   
    double sum1 = 0;
    double sum2 = 0;
  
    /* Calculation of L1 */
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (riskset(j,i) > 0) {
          sum1 = sum1 +
            delta(j,i)*I1(i,j)*deriv_vec[m](j,i)*(
                I3(i,j) -
                I2(i,j)*(std::exp(logtheta(j,i)))/(riskset(j,i) + I2(i,j)*(std::exp(logtheta(j,i)) - 1))
            );
        } else {sum1 = sum1 + 0;}
      }
    }
  
    /* Calculation of L2 */
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (riskset(i,j) > 0) {
          sum2 = sum2 +
            delta(i,j)*I2(j,i)*deriv_vec[m](i,j)*(
                I4(i,j) -
                I1(j,i)*(std::exp(logtheta(i,j)))/(riskset(i,j) + I1(j,i)*(std::exp(logtheta(i,j)) - 1))
            );
        } else {sum2 = sum2 + 0;}
      }
    }
  
    result[m] = -sum1-sum2;
    
  }

  return(result);
  
}


// [[Rcpp::export]]
NumericMatrix hessianC(const NumericMatrix riskset,
                const NumericMatrix logtheta,
                const Rcpp::List deriv,
                const int df,
                const NumericMatrix delta,
                const NumericMatrix I1,
                const NumericMatrix I2) {

  int n = riskset.nrow();
  int totalparam = df*df;

  NumericMatrix result(totalparam);

  /* Transform list of derivative matrices into vector of matrices */
  std::vector<NumericMatrix> deriv_vec(totalparam);
  for(int k = 0; k < totalparam; ++k) {
    NumericMatrix deriv_R = deriv[k];
    /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
    deriv_vec[k] = derivMat; */
    deriv_vec[k] = deriv_R;
  }

  for (int m = 0; m < totalparam; m++) {
    for (int l = m; l < totalparam; l++) {

      double sum1 = 0;
      double sum2 = 0;

      // L1
      for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
          if (riskset(j,i) > 0) {
            sum1 = sum1 -
              delta(j,i)*I1(i,j)*deriv_vec[m](j,i)*deriv_vec[l](j,i)*(riskset(j,i) - I2(i,j))*I2(i,j)*(std::exp(logtheta(j,i)))/pow(riskset(j,i) - I2(i,j) + I2(i,j)*(std::exp(logtheta(j,i))),2);
          } else {sum1 = sum1 + 0;}
        }
      }

      // L2
      for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
          if (riskset(i,j) > 0) {
            sum2 = sum2 -
              delta(i,j)*I2(j,i)*deriv_vec[m](i,j)*deriv_vec[l](i,j)*(riskset(i,j) - I1(j,i))*I1(j,i)*(std::exp(logtheta(i,j)))/pow(riskset(i,j) - I1(j,i) + I1(j,i)*(std::exp(logtheta(i,j))),2);
          } else {sum2 = sum2 + 0;}
        }
      }

      result(m,l) = -sum1-sum2;
      result(l,m) = result(m,l);

    }
  }

return(result);

}

// // [[Rcpp::export]]
// NumericVector MatToVec(const NumericMatrix matrix) {
// 
//   int row = matrix.nrow();
//   int column = matrix.ncol();
//   NumericVector vector(row*column);
// 
//   int k = 0;
// 
//   for (int j = 0; j < column; j++) {
//     for (int i = 0; i < row; i++) {
//       vector[k] = matrix(i,j);
//       k = k + 1;
//     }
//   }
// 
//   return(vector);
// 
// }
// 
// // [[Rcpp::export]]
// NumericMatrix VecToMat(const NumericVector vector) {
// 
//   int df = pow(vector.size(), 1/2);
//   NumericMatrix matrix(df);
// 
//   int k = 0;
//   for (int j = 0; j < df; j++) {
//     for (int i = 0; i < df; i++) {
//       matrix(i,j) = vector[k];
//       k = k + 1;
//     }
//   }
// 
//   return(matrix);
// 
// }



