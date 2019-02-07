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
    } else if(type==2){
      b=cest[12];
    } else{
      b=cest[14];
    }
    return b;
    
  } else if(let==2){
    
    if(type==1){
      b=cest[8];
    } else if(type==2){
      b=0;
    } else{
      b=0;
    }
    return b;
    
  } else if(let==3){
    
    if(type==1){
      b=cest[9];
    } else if(type==2){
      b=0;
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
      }else {
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
NumericVector CumH_c_ED(NumericVector v1, NumericVector v2, NumericVector u, NumericVector cest, 
                        NumericVector bref, bool affp, int type, NumericVector phi) {
  int q = u.size();
  NumericVector Hvec(2);
  
  NumericVector uv = u;
  double finalH1;
  double finalH0;
  double phisc=phi[0], phior=phi[1];
  
  
  for(int i=0; i < q; ++i) {
    
    double u = uv[i];
    
    NumericVector k1;
    k1 = v1[v1<u];
    int size = k1.size();
    
    if (size==0) {
      // double chaz(double lam1, double rho1, double betasc, double betaor, double phisc, double phior, 
      //             double lowsc, double lowor, double low, double upper)
      
      if(type==1){
        finalH0 = -chaz(cest[0],cest[1],0,0,0,0, 0,0, 0,u);
        finalH1 = finalH0;
      } else if (type==2){
        finalH0 = -chaz(cest[2],cest[3],0,0,0,0, 0,0, 0,u);
        finalH1 = finalH0;
      } else if (type==3){
        finalH0 = -chaz(cest[4],cest[5],0,0,0,0, 0,0, 0,u);
        finalH1 = finalH0;
      }
      
    } else {
      k1.push_back(u);
      size = k1.size();
      
      NumericVector intg, intw, H0full, H1full, v2his;
      double beta, beta_sc, beta_or, lowsc, lowor;
      beta_or = 0, beta_sc = 0, lowsc=0, lowor=0, beta=0;
      for (int i=0; i < size; ++i) {
        
        double beta; 
        if(i!=0){
          beta = coeff_c_ED(v2[i-1],type=type, affp, cest);
          if(v2[i-1]==1 || v2[i-1]==2 || v2[i-1]==3 ){
            beta_sc = beta;
            lowsc = k1[i-1];
          }else if(v2[i-1]==4){
            beta_or = beta;
            lowor = k1[i-1];
          }
        }
        
        if(beta==99){
          break;
        }
        
        // if(type!=1){
        //   
        //   double intg, intw = 0;
        //   
        // } else {
        //   
        //   NumericVector inter(2);
        //   
        //   if(i==0){
        //     inter[0]=0;
        //     inter[1]=0;
        //     
        //   } else {
        //     inter = coeff_interaction_c_IS(v2[i-1], v2his, bref);
        //     v2his.push_back(v2[i-1]);  
        //   }
        //   
        //   intw.push_back(inter[1]);
        //   intg.push_back(inter[0]);
        // }
        //
        double Hpiece, H_carrier, H_noncarrier, low, upper;
        
        upper = k1[i];
        if(i!=0){
        low = k1[i-1];   
        } else{
        low = 0;
        }
        // double chaz(double lam1, double rho1, double betasc, double betaor, double phisc, double phior, 
        //             double lowsc, double lowor, double low, double upper)
          
        if(type==1){
         Hpiece = chaz(cest[0],cest[1],beta_sc,beta_or,phisc,phior,lowsc,lowor, low,upper);
        } else if (type==2){
         Hpiece = chaz(cest[2],cest[3],beta_sc,0,0,0, 0,0, low,upper);
        } else {
         Hpiece = chaz(cest[4],cest[5],beta_sc,0,0,0, 0,0, low,upper);
        }
        
        H_carrier = Hpiece;
        H_noncarrier = Hpiece;
        
        H1full.push_back(H_carrier);
        H0full.push_back(H_noncarrier);
      }
      finalH1 = -sum(H1full);
      finalH0 = -sum(H0full);
    }
    Hvec[0] = finalH0;
    Hvec[1] = finalH1;
  }
  return Hvec;
}

/*** R
microbenchmark(CumH_c_IS(v1=v1,v2=v2,u=time0,cest=cest,bref=bref,affp=FALSE,type=1),times=1)
microbenchmark(CumH_c_ED(v1=v1,v2=v2,u=time0,cest=cest,bref=bref,affp=FALSE,type=1,phi=phiv), times=1)
*/
