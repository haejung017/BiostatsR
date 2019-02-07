penf <- function(est, age, sex, mut, base.dist="Weibull", frailty.dist=NULL, agemin){
  if(base.dist=="lognormal") base.est <- c(est[1], exp(est[2]))
  else base.est <- exp(est[1:2])
  H <- cumhaz(base.dist, age-agemin, base.est)*exp(est[3]*sex+est[4]*mut) 
  if(is.null(frailty.dist)) pen = 1-exp(-H)
  else {
    if(length(est)!=5) stop("length(est) should be 5.")
    else pen = 1-laplace(dist=frailty.dist, g=H, k=est[5])
  }
  names(pen) <- NULL
  return(pen)
}