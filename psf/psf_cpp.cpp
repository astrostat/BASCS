#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector psf_cpp(double offangle, double ellip, double slope, double psfnorm, double ro, NumericVector pts, NumericVector center)
{
  int m = pts.size()/2;
  double x[m]; double y[m]; double ci; double si; double rad[m];
  for (int i=0; i<m; i++){
    x[i] = pts[i] - center[0];
    y[i] = pts[i+m] - center[1];
  }
  ci = cos(offangle);
  si = sin(offangle);
  for (int i=0; i<m; i++){
    rad[i] = sqrt (pow(x[i]*ci+y[i]*si,2) + pow(y[i]*ci-x[i]*si,2)/pow(1-ellip,2));
  }
  Rcpp::NumericVector val(m);
  for (int i=0; i<m; i++){
    val[i] = psfnorm/pow(1+pow(rad[i]/ro,2),slope);
  }
  return val;
}

