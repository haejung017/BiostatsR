## under Permanent Exposure (PE) TVC
cH_PE<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  
  bH0 <- function(u){  # Cumulative hazard function
    (lamb*u)^rho
  }
  
  
  return( exp(bz)*(bH0(b) - bH0 (a) ) )
  
}

## under Exponential Decay (ED) TVC
cH_ED<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  phi<-parm[4]
  
  bh0<-function(u){
    (lamb*rho)*(lamb*u)^(rho-1)
  }
  
  integrand <- function(u) { bh0(u)*exp(bz*exp(-(u-a)*phi)) }
  int <- try( integrate( integrand, lower=a, upper=b), silent = TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <-NA_real_
  } else {
    integrated <-int$value
  }
  
  return(integrated)
  
}


## under Cox and Oaks (CO) TVC
cH_CO<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  phi<-parm[4]
  a0<-parm[5]
  
  bh0<-function(u){
    (lamb*rho)*(lamb*u)^(rho-1)
  }
  
  integrand <- function(u) { bh0( u)*exp(a0+bz*exp(-(u-a)*phi)) }
  int <- try( integrate( integrand, lower=a, upper=b), silent = TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <-NA_real_
  } else {
    integrated <-int$value
  }
  
  return(integrated)
  
}

# combine all these function in one for subject j
cHj<-function(x, parm,TDfunc){
  a=x[1]; b=x[2]
  if(TDfunc=="PE") res<-cH_PE(a, b, parm)
  else if(TDfunc=="ED") res<-cH_ED(a, b, parm)
  else if(TDfunc=="CO") res<-cH_CO(a, b, parm)
  else stop("TDfunc should be PE or ED or CO")
  res
}

#for several subjects :  use lapply function
cH<-function(A, B, parm, TDfunc){
  unlist( lapply(1:length(A), function(i) cHj(x=c(A[i],B[i]), parm=parm, TDfunc=TDfunc) ) )
}
