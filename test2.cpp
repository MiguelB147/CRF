#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix IndGreater(NumericVector &x) {
  int n = x.size();
  NumericMatrix elem(n);
  for (int j=0; j<n; j++) {
    for (int i=0; i<n; i++) {
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
NumericMatrix IndLess(NumericVector &x) {
  int n = x.size();
  NumericMatrix elem(n);
  for (int j=0; j<n; j++) {
    for (int i=0; i<n; i++) {
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
NumericMatrix IndEqual(NumericVector &x) {
  int n = x.size();
  NumericMatrix elem(n);
  for (int j=0; j<n; j++) {
    for (int i=0; i<n; i++) {
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
int Ind2(NumericVector &x, NumericVector &y, double &a, double &b) {
  int n = x.size();
  int sum = 0;
  for (int i=0; i<n; i++) {
    if (x[i] >= a and y[i] >= b) {
      sum += 1;
    } else {
      sum += 0;
    }
  }
  return sum;
}

// [[Rcpp::export]]
NumericMatrix risksetC(NumericVector &x, NumericVector &y) {
  int n = x.size();
  NumericMatrix riskset(n);
  for (int j=0; j<n; j++) {
    for (int i=0; i<n; i++) {
      riskset (j,i) = Ind2(x, y, x[j], y[i]);
    }
  }
  return riskset;
}

// [[Rcpp::export]]
NumericMatrix DeltaC(NumericVector &x, NumericVector &y) {
  int n = x.size();
  NumericMatrix delta(n);
  for (int j=0; j<n; j++) {
    for (int i=0; i<n; i++) {
      delta (j,i) = x[j]*y[i];
    }
  }
  return delta;
}



// [[Rcpp::export]]
double logLikC(const NumericVector &riskset1,
               const NumericVector &riskset2,
               const NumericVector &logtheta1,
               const NumericVector &logtheta2,
               const NumericVector &delta1,
               const NumericVector &delta2,
               const NumericVector &I1,
               const NumericVector &I2,
               const NumericVector &I3,
               const NumericVector &I4,
               const NumericVector &I5,
               const NumericVector &I6) {

  double sum1;
  double sum2;
  
  sum1 = Rcpp::sum(delta1*I1*(logtheta1*I5 - Rcpp::log(riskset1 + I2*Rcpp::exp(logtheta1) - I2)));
  sum2 = Rcpp::sum(delta2*I3*(logtheta2*I6 - Rcpp::log(riskset2 + I4*Rcpp::exp(logtheta2) - I4)));
  
  return(-sum1-sum2);
}

// [[Rcpp::export]]
NumericVector gradientC(const NumericVector &riskset1,
                        const NumericVector &riskset2,
                        const NumericVector &logtheta1,
                        const NumericVector &logtheta2,
                        const Rcpp::List &deriv,
                        const int &df,
                        const NumericVector &delta1,
                        const NumericVector &delta2,
                        const NumericVector &I1,
                        const NumericVector &I2,
                        const NumericVector &I3,
                        const NumericVector &I4,
                        const NumericVector &I5,
                        const NumericVector &I6
) {

  int n = riskset1.length();
  int totalparam = df*df;
  NumericVector result(totalparam);
  NumericVector common1(n);
  NumericVector common2(n);

  /* Transform list of derivative matrices into vector of matrices */
  std::vector<NumericMatrix> deriv_vec(totalparam);
  for (int k = 0; k < totalparam; ++k) {
    NumericMatrix deriv_R = deriv[k];
    deriv_vec[k] = deriv_R;
  }
  

  common1 = delta1*I1*(I5 - I2*Rcpp::exp(logtheta1)/(riskset1 + I2*Rcpp::exp(logtheta1) - I2));
  common2 = delta2*I3*(I6 - I4*Rcpp::exp(logtheta2)/(riskset2 + I4*Rcpp::exp(logtheta2) - I4));

  for (int m=0; m<totalparam; m++) {
    
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    /* Calculation of L1 */
    for (int j=0; j<n; j++) {
      sum1 += common1(j)*deriv_vec[m](j,0);
      sum2 += common2(j)*deriv_vec[m](j,1);
    }
    
    result(m) = -sum1-sum2;
    
  }

return(result);

}


// [[Rcpp::export]]
NumericVector gradientPoly(const NumericVector &riskset1,
                           const NumericVector &riskset2,
                           const NumericVector &logtheta1,
                           const NumericVector &logtheta2,
                           const Rcpp::List &deriv,
                           const int &df,
                           const NumericVector &delta1,
                           const NumericVector &delta2,
                           const NumericVector &I1,
                           const NumericVector &I2,
                           const NumericVector &I3,
                           const NumericVector &I4,
                           const NumericVector &I5,
                           const NumericVector &I6) {
  
  int n = riskset1.length();
  NumericVector result(df);
  NumericVector common1(n);
  NumericVector common2(n);
  
  /* Transform list of derivative matrices into vector of matrices */
  std::vector<NumericMatrix> deriv_vec(df);
  for(int k = 0; k < df; ++k) {
    NumericMatrix deriv_R = deriv[k];
    /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
     deriv_vec[k] = derivMat; */
    deriv_vec[k] = deriv_R;
  }

  common1 = delta1*I1*(I5 - I2*Rcpp::exp(logtheta1)/(riskset1 + I2*Rcpp::exp(logtheta1) - I2));
  common2 = delta2*I3*(I6 - I4*Rcpp::exp(logtheta2)/(riskset2 + I4*Rcpp::exp(logtheta2) - I4));
  
  for (int m=0; m<df; m++) {
    
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    /* Calculation of L1 */
    for (int j=0; j<n; j++) {
      sum1 += common1(j)*deriv_vec[m](j,0);
      sum2 += common2(j)*deriv_vec[m](j,1);
    }
  
    result(m) = -sum1-sum2;
    
  }

  return(result);
  
}


// [[Rcpp::export]]
NumericMatrix hessianC(const NumericVector &riskset1,
                       const NumericVector &riskset2,
                       const NumericVector &logtheta1,
                       const NumericVector &logtheta2,
                       const Rcpp::List &deriv,
                       const int &df,
                       const NumericVector &delta1,
                       const NumericVector &delta2,
                       const NumericVector &I1,
                       const NumericVector &I2,
                       const NumericVector &I3,
                       const NumericVector &I4) {

  int n = riskset1.length();
  int totalparam = df*df;
  NumericVector common1(n);
  NumericVector common2(n);

  NumericMatrix result(totalparam);

  /* Transform list of derivative matrices into vector of matrices */
  std::vector<NumericMatrix> deriv_vec(totalparam);
  for(int k = 0; k < totalparam; ++k) {
    NumericMatrix deriv_R = deriv[k];
    /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
    deriv_vec[k] = derivMat; */
    deriv_vec[k] = deriv_R;
  }
  
    common1 = -delta1*I1*(riskset1 - I2)*I2*Rcpp::exp(logtheta1)/Rcpp::pow(riskset1 - I2 + I2*Rcpp::exp(logtheta1),2);
    common2 = -delta2*I3*(riskset2 - I4)*I4*Rcpp::exp(logtheta2)/Rcpp::pow(riskset2 - I4 + I4*Rcpp::exp(logtheta2),2);  

  for (int m = 0; m < totalparam; m++) {
    for (int l = m; l < totalparam; l++) {
      
      double sum1 = 0.0;
      double sum2 = 0.0;

      for (int j=0; j<n; j++) {
        sum1 += common1(j)*deriv_vec[m](j,0)*deriv_vec[l](j,0);
        sum2 += common2(j)*deriv_vec[m](j,1)*deriv_vec[l](j,1);
      }

      result(m,l) = -sum1-sum2;
      result(l,m) = result(m,l);

    }
  }

return(result);

}

// [[Rcpp::export]]
NumericMatrix hessianPolyC(const NumericVector &riskset1,
                           const NumericVector &riskset2,
                           const NumericVector &logtheta1,
                           const NumericVector &logtheta2,
                           const Rcpp::List &deriv,
                           const int &df,
                           const NumericVector &delta1,
                           const NumericVector &delta2,
                           const NumericVector &I1,
                           const NumericVector &I2,
                           const NumericVector &I3,
                           const NumericVector &I4) {
  
  int n = riskset1.length();
  NumericVector common1(n);
  NumericVector common2(n);

  NumericMatrix result(df);

  /* Transform list of derivative matrices into vector of matrices */
  std::vector<NumericMatrix> deriv_vec(df);
  for(int k = 0; k < df; ++k) {
    NumericMatrix deriv_R = deriv[k];
    /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
    deriv_vec[k] = derivMat; */
    deriv_vec[k] = deriv_R;
  }
  
  common1 = -delta1*I1*(riskset1 - I2)*I2*Rcpp::exp(logtheta1)/Rcpp::pow(riskset1 - I2 + I2*Rcpp::exp(logtheta1),2);
  common2 = -delta2*I3*(riskset2 - I4)*I4*Rcpp::exp(logtheta2)/Rcpp::pow(riskset2 - I4 + I4*Rcpp::exp(logtheta2),2);  
  
  for (int m = 0; m < df; m++) {
    for (int l = m; l < df; l++) {
      
      double sum1 = 0.0;
      double sum2 = 0.0;
    
      for (int j=0; j<n; j++) {
        sum1 += common1(j)*deriv_vec[m](j,0)*deriv_vec[l](j,0);
        sum2 += common2(j)*deriv_vec[m](j,1)*deriv_vec[l](j,1);
      }
    
      result(m,l) = -sum1-sum2;
      result(l,m) = result(m,l);
    
    }
  }
  
  return(result);
  
}

NumericMatrix testfunct(const List &deriv,
                        const int &df) {
  
  int totalparam = df*df;
  
  std::vector<NumericMatrix> deriv_vec(totalparam);
  for(int k = 0; k < totalparam; ++k) {
    NumericMatrix deriv_R = deriv[k];
    /* arma::mat derivMat(deriv_R.begin(), deriv_R.rows(), deriv_R.cols(), false, true);
    deriv_vec[k] = derivMat; */
    deriv_vec[k] = deriv_R;
  }
  NumericMatrix result = deriv_vec[0];
  return(result);
}

