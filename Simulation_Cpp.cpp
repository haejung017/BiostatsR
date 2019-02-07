// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include <Rcpp.h>
namespace Numer {}
using namespace Numer;
using namespace Rcpp;

class hazard: public Func
{
private:
  double rho1;
  double betasc;
  double phisc;
  double lowsc;
public:
  hazard(double rho1_, double betasc_, double phisc_, double lowsc_) : rho1(rho1_),
  betasc(betasc_), phisc(phisc_), lowsc(lowsc_) {}
  
  double operator()(const double& x) const
  {
    return pow(x,(rho1-1))*exp(betasc*exp(-(x-lowsc)*phisc));
  }
};


// [[Rcpp::export]]
double cumulhazC(double rho1, double betasc, double phisc, 
            double lowsc, double lower, double upper)
{
  // const double true_val = R::pbeta(upper, a, b, 1, 0) -
  //   R::pbeta(lower, a, b, 1, 0);
  
  hazard f(rho1, betasc, phisc, lowsc);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return res;
}
