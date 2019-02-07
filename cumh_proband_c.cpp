#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]
bool contains(NumericVector X, double z) {
  return std::find(X.begin(), X.end(), z)!=X.end();
}
// [[Rcpp::export]]
double H1_c(double t, NumericVector est) {
  
  double cumh;
  double lambda1 = est[0];
  double rho1 = est[1];
  cumh = pow(lambda1*t,rho1);
  
  return cumh; 
  
}
// [[Rcpp::export]]
double H2_c(double t, NumericVector est) {
  
  double cumh;
  double lambda2 = est[2];
  double rho2 = est[3];
  cumh = pow(lambda2*t,rho2);
  
  return cumh; 
  
}
// [[Rcpp::export]]
double H3_c(double t, NumericVector est) {
  
  double cumh;
  double lambda3 = est[4];
  double rho3 = est[5];
  cumh = pow(lambda3*t,rho3);
  
  return cumh; 
  
}
// [[Rcpp::export]]
double coeff_c(int let, int type, bool affp, NumericVector cest) {
  
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

// [[Rcpp::export]]
NumericVector coeff_interaction_c(int let, NumericVector v2his, NumericVector bref) {

    NumericVector out(2);
    
    if(let==5){
      out[0] = 0;
      out[1] = 0;
    } else if(let==4){
      if(contains(v2his,1)){
        out[0] =0;
        out[1] =bref[3];
        if(contains(v2his,2)){
          out[0]=0;
          out[1]= bref[3]+bref[4];
          if(contains(v2his,3)){
            out[0] =0;
            out[1] = bref[3]+bref[4]+bref[5];
          }
        }
      } else{
        out[0] =0;
        out[1] =0;
      }
    } else{
      if(contains(v2his,4)){
        out[0] = bref[let-1];
        out[1] = bref[let+2];
      } else{
        out[0] = bref[let-1];
        out[1] = 0;
      }
    }
    return out;
  }
  
// [[Rcpp::export]]
NumericVector CumH_proband_c(NumericVector v1, NumericVector v2, NumericVector u, NumericVector cest, NumericVector bref, bool affp, int type) {
    int q = u.size();
    NumericVector Hvec(q);

    NumericVector uv = u;
    double finalH;


    for(int i=0; i < q; ++i) {

      double u = uv[i];

      NumericVector k1;
      k1 = v1[v1<u];
      int size = k1.size();

       if (size==0) {

        if(type==1){
          finalH = -H1_c(u,cest);
        } else if (type==2){
          finalH = -H2_c(u,cest);
        } else if (type==3){
          finalH = -H3_c(u,cest);
        }

       } else {
            k1.push_back(u);
            size = k1.size();

            NumericVector intg, intw, H0full, H1full, beta_vec, v2his;

                for (int i=0; i < size; ++i) {

                double beta; 
                if(i==0){
                  beta= 0;
                }else{
                  // int v2e;
                  // v2e = v2[i-1];
                  beta = coeff_c(v2[i-1],type=type, affp, cest);
                }
                
                 if(beta==99){
                  break;
                 }

                beta_vec.push_back(beta);

                 if(type!=1){

                  double intg, intw = 0;
                   
                 } else {
                  
                  NumericVector inter(2);
                  
                  if(i==0){
                  inter[0]=0;
                  inter[1]=0;
                    
                  } else {
                  inter = coeff_interaction_c(v2[i-1], v2his, bref);
                  v2his.push_back(v2[i-1]);  
                  }
                  
                  intw.push_back(inter[1]);
                  intg.push_back(inter[0]);
                }
          //
                double Hpiece_u, Hpiece_l, H_carrier;
                
                double up = k1[i];
                
                if(type==1){
                Hpiece_u = H1_c(up,cest);
                  if(i==0){
                  Hpiece_l = 0;
                  } else{
                  double low = k1[i-1];  
                  Hpiece_l = H1_c(low,cest);  
                  }
                } else if (type==2){
                Hpiece_u = H2_c(up,cest);
                  if(i==0){
                    Hpiece_l = 0;
                  } else{
                    double low = k1[i-1];
                    Hpiece_l = H2_c(low,cest);  
                  }
                } else if (type==3){
                Hpiece_u = H3_c(up,cest);
                  if(i==0){
                    Hpiece_l = 0;
                  } else{
                    double low = k1[i-1];
                    Hpiece_l = H3_c(low,cest);  
                  }
                }

                H_carrier = (Hpiece_u - Hpiece_l)*exp(sum(beta_vec)+sum(intw)+sum(intg));

                H1full.push_back(H_carrier);
                }
              finalH = -sum(H1full);
              }
           Hvec[i] = finalH;
           }
   return Hvec;
}

// [[Rcpp::export]]
NumericVector CumH_c(NumericVector v1, NumericVector v2, NumericVector u, NumericVector cest, NumericVector bref, bool affp, int type) {
  int q = u.size();
  NumericVector Hvec(2);
  
  NumericVector uv = u;
  double finalH1;
  double finalH0;
  
  
  for(int i=0; i < q; ++i) {
    
    double u = uv[i];
    
    NumericVector k1;
    k1 = v1[v1<u];
    int size = k1.size();
    
    if (size==0) {
      
      if(type==1){
        finalH0 = -H1_c(u,cest);
        finalH1 = -H1_c(u,cest);
      } else if (type==2){
        finalH0 = -H2_c(u,cest);
        finalH1 = -H2_c(u,cest);
      } else if (type==3){
        finalH0 = -H3_c(u,cest);
        finalH1 = -H3_c(u,cest);
      }
      
    } else {
      k1.push_back(u);
      size = k1.size();
      
      NumericVector intg, intw, H0full, H1full, beta_vec, v2his;
      
      for (int i=0; i < size; ++i) {
        
        double beta; 
        if(i==0){
          beta= 0;
        }else{
          // int v2e;
          // v2e = v2[i-1];
          beta = coeff_c(v2[i-1],type=type, affp, cest);
        }
        
        if(beta==99){
          break;
        }
        
        beta_vec.push_back(beta);
        
        if(type!=1){
          
          double intg, intw = 0;
          
        } else {
          
          NumericVector inter(2);
          
          if(i==0){
            inter[0]=0;
            inter[1]=0;
            
          } else {
            inter = coeff_interaction_c(v2[i-1], v2his, bref);
            v2his.push_back(v2[i-1]);  
          }
  
          intw.push_back(inter[1]);
          intg.push_back(inter[0]);
        }
        //
        double Hpiece_u, Hpiece_l, H_carrier, H_noncarrier;
        
        double up = k1[i];
        
        if(type==1){
          Hpiece_u = H1_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];  
            Hpiece_l = H1_c(low,cest);  
          }
        } else if (type==2){
          Hpiece_u = H2_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];
            Hpiece_l = H2_c(low,cest);  
          }
        } else if (type==3){
          Hpiece_u = H3_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];
            Hpiece_l = H3_c(low,cest);  
          }
        }
        
        H_carrier = (Hpiece_u - Hpiece_l)*exp(sum(beta_vec)+sum(intw)+sum(intg));
        H_noncarrier = (Hpiece_u - Hpiece_l)*exp(sum(beta_vec)+sum(intw));
        
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


// Interval specific screening effects 

// [[Rcpp::export]]
double coeff_c_IS(int let, int type, bool affp, NumericVector cest) {
  
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
      b=cest[8]-cest[7];
    } else if(type==2){
      b=0;
    } else{
      b=0;
    }
    return b;
    
  } else if(let==3){
    
    if(type==1){
      b=cest[9]-cest[8];
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

// [[Rcpp::export]]
NumericVector coeff_interaction_c_IS(int let, NumericVector v2his, NumericVector bref) {
  
  NumericVector out(2);
  
  if(let==5){
    out[0] = 0;
    out[1] = 0;
  } else if(let==4){ // encounter oophorectomy, check retrospectively screening activities
    if(contains(v2his,1)){
      out[0] =0;
      out[1] =bref[3];
      if(contains(v2his,2)){
        out[0]=0;
        out[1]=bref[4];
        if(contains(v2his,3)){
          out[0] =0;
          out[1] =bref[5];
        }
      }
    } else{
      out[0] =0;
      out[1] =0;
    }
  } else if(let==1){ // encounter screenings
    if(contains(v2his,4)){
      out[0] = bref[0];
      out[1] = bref[3];
      } else{
      out[0] = bref[0];
      out[1] = 0;
     }
  } else if(let==2){
    if(contains(v2his,4)){
      out[0] = bref[1]-bref[0];
      out[1] = bref[4]-bref[3];
    } else{
      out[0] = bref[1]-bref[0];
      out[1] = 0;
    }    
  } else if(let==3){
    if(contains(v2his,4)){
      out[0] = bref[2]-bref[1];
      out[1] = bref[5]-bref[4];
    } else{
      out[0] = bref[2]-bref[1];
      out[1] = 0;
    }    
  }
  return out;
}

// [[Rcpp::export]]
NumericVector CumH_c_IS(NumericVector v1, NumericVector v2, NumericVector u, NumericVector cest, NumericVector bref, bool affp, int type) {
  int q = u.size();
  NumericVector Hvec(2);
  
  NumericVector uv = u;
  double finalH1;
  double finalH0;
  
  
  for(int i=0; i < q; ++i) {
    
    double u = uv[i];
    
    NumericVector k1;
    k1 = v1[v1<u];
    int size = k1.size();
    
    if (size==0) {
      
      if(type==1){
        finalH0 = -H1_c(u,cest);
        finalH1 = -H1_c(u,cest);
      } else if (type==2){
        finalH0 = -H2_c(u,cest);
        finalH1 = -H2_c(u,cest);
      } else if (type==3){
        finalH0 = -H3_c(u,cest);
        finalH1 = -H3_c(u,cest);
      }
      
    } else {
      k1.push_back(u);
      size = k1.size();
      
      NumericVector intg, intw, H0full, H1full, beta_vec, v2his;
      
      for (int i=0; i < size; ++i) {
        
        double beta; 
        if(i==0){
          beta= 0;
        }else{
          beta = coeff_c_IS(v2[i-1],type=type, affp, cest);
        }
        
        if(beta==99){
          break;
        }
        
        beta_vec.push_back(beta);
        
        if(type!=1){
          
          double intg, intw = 0;
          
        } else {
          
          NumericVector inter(2);
          
          if(i==0){
            inter[0]=0;
            inter[1]=0;
            
          } else {
            inter = coeff_interaction_c_IS(v2[i-1], v2his, bref);
            v2his.push_back(v2[i-1]);  
          }
          
          intw.push_back(inter[1]);
          intg.push_back(inter[0]);
        }
        //
        double Hpiece_u, Hpiece_l, H_carrier, H_noncarrier;
        
        double up = k1[i];
        
        if(type==1){
          Hpiece_u = H1_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];  
            Hpiece_l = H1_c(low,cest);  
          }
        } else if (type==2){
          Hpiece_u = H2_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];
            Hpiece_l = H2_c(low,cest);  
          }
        } else if (type==3){
          Hpiece_u = H3_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];
            Hpiece_l = H3_c(low,cest);  
          }
        }
        
        H_carrier = (Hpiece_u - Hpiece_l)*exp(sum(beta_vec)+sum(intw)+sum(intg));
        H_noncarrier = (Hpiece_u - Hpiece_l)*exp(sum(beta_vec)+sum(intw));
        
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

// [[Rcpp::export]]
NumericVector CumH_proband_c_IS(NumericVector v1, NumericVector v2, NumericVector u, NumericVector cest, NumericVector bref, bool affp, int type) {
  int q = u.size();
  NumericVector Hvec(q);
  
  NumericVector uv = u;
  double finalH;
  
  
  for(int i=0; i < q; ++i) {
    
    double u = uv[i];
    
    NumericVector k1;
    k1 = v1[v1<u];
    int size = k1.size();
    
    if (size==0) {
      
      if(type==1){
        finalH = -H1_c(u,cest);
      } else if (type==2){
        finalH = -H2_c(u,cest);
      } else if (type==3){
        finalH = -H3_c(u,cest);
      }
      
    } else {
      k1.push_back(u);
      size = k1.size();
      
      NumericVector intg, intw, H0full, H1full, beta_vec, v2his;
      
      for (int i=0; i < size; ++i) {
        
        double beta; 
        if(i==0){
          beta= 0;
        }else{
          // int v2e;
          // v2e = v2[i-1];
          beta = coeff_c_IS(v2[i-1],type=type, affp, cest);
        }
        
        if(beta==99){
          break;
        }
        
        beta_vec.push_back(beta);
        
        if(type!=1){
          
          double intg, intw = 0;
          
        } else {
          
          NumericVector inter(2);
          
          if(i==0){
            inter[0]=0;
            inter[1]=0;
            
          } else {
            inter = coeff_interaction_c_IS(v2[i-1], v2his, bref);
            v2his.push_back(v2[i-1]);  
          }
          
          intw.push_back(inter[1]);
          intg.push_back(inter[0]);
        }
        //
        double Hpiece_u, Hpiece_l, H_carrier;
        
        double up = k1[i];
        
        if(type==1){
          Hpiece_u = H1_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];  
            Hpiece_l = H1_c(low,cest);  
          }
        } else if (type==2){
          Hpiece_u = H2_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];
            Hpiece_l = H2_c(low,cest);  
          }
        } else if (type==3){
          Hpiece_u = H3_c(up,cest);
          if(i==0){
            Hpiece_l = 0;
          } else{
            double low = k1[i-1];
            Hpiece_l = H3_c(low,cest);  
          }
        }
        
        H_carrier = (Hpiece_u - Hpiece_l)*exp(sum(beta_vec)+sum(intw)+sum(intg));
        
        H1full.push_back(H_carrier);
      }
      finalH = -sum(H1full);
    }
    Hvec[i] = finalH;
  }
  return Hvec;
}

