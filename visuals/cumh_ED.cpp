// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include <Rcpp.h>
namespace Numer {}
using namespace Numer;
using namespace Rcpp;

// P(0.3 < X < 0.8), X ~ Beta(a, b)
// [[Rcpp::export]]
double coeff_c_ED(int let, int type, bool affp, NumericVector cest) {
  
  double b;
  
  if(let==1){
    if(type==1){
      b=cest[7];
    } else {
      b=0;
    }
    return b;
    
  } else if(let==2){
    
    if(type==1){
      b=cest[8];
    } else{
      b=0;
    }
    return b;
    
  } else if(let==3){
    
    if(type==1){
      b=cest[9];
    } else{
      b=0;
    }
    return b;
    
  } else if(let==4){
    
    if(type==1){
      b=cest[10];
    } else if(type==2){
      if(affp==FALSE){
        b=99;
      }else{
        b=0;    
      }
    } else{
      b=0;
    }
    return b;
    
  } else {
    if(type==1){
      if(affp==FALSE){
        b = 99;
      }else{
        b = 0;
      }
    }else{
      b=0;
    }
    return b;
  }
}

//(lambda1^rho1)*rho1*time0^(rho1-1)
//
class hazard: public Func
{
private:
  double lam1;
  double rho1;
  double betasc;
  double betaor;
  double phisc;
  double phior;
  double lowsc;
  double lowor;
public:
  hazard(double lam1_, double rho1_, double betasc_, double betaor_, 
         double phisc_, double phior_, double lowsc_,double lowor_) : lam1(lam1_), rho1(rho1_),
         betasc(betasc_), betaor(betaor_), phisc(phisc_), phior(phior_), lowsc(lowsc_), lowor(lowor_) {}
  
  double operator()(const double& x) const
  {
    return pow(lam1,rho1)*rho1*pow(x,(rho1-1))*exp(betasc*exp(-(x-lowsc)*phisc) + betaor*exp(-(x-lowor)*phior));
  }
};


// [[Rcpp::export]]
double chaz(double lam1, double rho1, double betasc, double betaor, double phisc, double phior, 
            double lowsc, double lowor, double low, double upper)
{
  const double lower = low;
  // const double true_val = R::pbeta(upper, a, b, 1, 0) -
  //   R::pbeta(lower, a, b, 1, 0);
  
  hazard f(lam1, rho1, betasc, betaor, phisc, phior, lowsc, lowor);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return res;
}

// [[Rcpp::export]]
NumericMatrix CumH_c_ED(NumericMatrix v1, NumericMatrix v2, NumericVector u, NumericVector cest, 
                        bool affp, int type, NumericVector phi) {
  int n = v1.nrow();
  NumericMatrix Hvec(n,1);
  NumericVector uv = u;
  
  double finalH=0;
  double phisc=0, phior=phi[3];
  
  for(int j=0; j < n; ++j) {
    
    double u = uv[j];
    NumericVector k1 = v1( j , _ );
    NumericVector k2 = k1[k1<u];
    int size = k2.size();
    
    if (size==0) {
      // double chaz(double lam1, double rho1, double betasc, double betaor, double phisc, double phior, 
      //             double lowsc, double lowor, double low, double upper)
      
      if(type==1){
        finalH = -chaz(cest[0],cest[1],0,0,0,0, 0,0, 0,u);
      } else if (type==2){
        finalH = -chaz(cest[2],cest[3],0,0,0,0, 0,0, 0,u);
      } else if (type==3){
        finalH = -chaz(cest[4],cest[5],0,0,0,0, 0,0, 0,u);
      }
      
    } else {
      k2.push_back(u);
      size = k2.size();
      NumericVector w2 = v2( j , _ );
      
      NumericVector Hfull;
      double beta_sc=0, beta_or=0, lowsc=0, lowor=0;
      for (int i=0; i < size; ++i) {
        
        double beta; 
        if(i!=0){
          beta = coeff_c_ED(w2[i-1],type=type, affp, cest);
          if(w2[i-1]==1 || w2[i-1]==2 || w2[i-1]==3 ){
            beta_sc = beta;
            phisc = phi[w2[i-1]-1]; 
            lowsc = k2[i-1];
          }else if(w2[i-1]==4){
            beta_or = beta;
            lowor = k2[i-1];
          }
        }else{
          beta = 0;
          beta_sc = beta;
          beta_or = beta;
        }
        
        if(beta==99){
          break;
        }

        double Hpiece, H_noncarrier, low, upper;
        
        upper = k2[i];
        if(i!=0){
          low = k2[i-1];   
        } else{
          low = 0;
        }
        // double chaz(double lam1, double rho1, double betasc, double betaor, double phisc, double phior, 
        //             double lowsc, double lowor, double low, double upper)
        
        if(type==1){
          Hpiece = chaz(cest[0],cest[1],beta_sc,beta_or,phisc,phior,lowsc=lowsc,lowor=lowor, low=low,upper=upper);
        } else if (type==2){
          Hpiece = chaz(cest[2],cest[3],beta_sc,0,0,0, 0,0, low=low,upper=upper);
        } else {
          Hpiece = chaz(cest[4],cest[5],beta_sc,0,0,0, 0,0, low=low,upper=upper);
        }
        
        H_noncarrier = Hpiece;
        
        Hfull.push_back(H_noncarrier);
      }
      finalH = -sum(Hfull);
    }
    Hvec(j,0) = finalH;
  }
  return Hvec;
}
