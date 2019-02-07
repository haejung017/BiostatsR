summary.penmodel <- function(object, correlation=FALSE, ...){
   
  ans <- object[c("coefficients")]
  class(ans) <- "summary.penmodel"
  tval <-  object$coef/object$se
  pval <- 2 * pt(abs(tval), 1, lower.tail = FALSE)
  ans$coefficients <- cbind(object$coef, object$se, tval, pval)
  dimnames(ans$coefficients) <- list(names(object$coef), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  
  ans$pen70 <- rbind(Estimate=object$pen70.est, SE=object$pen70.se, CI=object$pen70.ci)
  ans$varcov <- object$varcov
  if(correlation) ans$correlation <- object$varcov/object$se^2
  return(ans) 


}
