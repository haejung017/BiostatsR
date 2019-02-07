
plot.penmodel <- function(x, print=TRUE, mark.time=FALSE, conf.int=FALSE, ...){

  base.dist <- attr(x,"base.dist")
  agemin <- attr(x, "agemin")
  if(base.dist=="lognormal") parms <- c(x$coefficients[1], exp(x$coefficients[2]))
  else parms <- exp(x$coefficients[1:2])
 penout <- penplot(base.parms=parms, vbeta=x$coefficients[3:4], base.dist=base.dist, variation="none", agemin=agemin, print=FALSE, ...)

 penest<-rbind(penout$pen[[1]],penout$pen[[2]],penout$pen[[2]],penout$pen[[4]])
 row.names(penest) <- names(x$pen70.est)
 data <- attr(x, "data")

 sfit<-survfit(Surv(time, status)~gender+mgene, data=data[data$proband==0 & !is.na(data$mgene),])
 lines(sfit, fun="event", lty=c(2,1,2,1), col=c("red","red","blue","blue"),
       mark.time=mark.time, conf.int=FALSE, ...)
 if(conf.int){
   lines(sfit, fun="event",
         conf.int="only", col=c("red","red","blue","blue"), lty=c(4,3,4,3))
   cat("Calculating ... \n")
   xx = penout$x.age
   ci <- sapply(xx, pen.ci, est=x$coefficients, cov=x$varcov, base.dist=base.dist, agemin=agemin, n=1000)
     lower<-ci[seq(3,16,4),]/100
     upper<-ci[seq(4,16,4),]/100
     row.names(lower)<-row.names(upper) <- names(x$pen70.est)
     lines(xx, lower[1,], lty=3, col="blue")
     lines(xx, upper[1,], lty=3, col="blue")
     lines(xx, lower[2,], lty=3, col="red")
     lines(xx, upper[2,], lty=3, col="red")
     lines(xx, lower[3,], lty=4, col="blue")
     lines(xx, upper[3,], lty=4, col="blue")
     lines(xx, lower[4,], lty=4, col="red")
     lines(xx, upper[4,], lty=4, col="red")
     
     legend("topright", c("95% CI","95% CI", "95% CI","95% CI"),
           col=c("blue","red"), lty=c(3,3,4,4), bty="n")

     out <- list(coefficients=x$coefficients, pen70=x$pen70.est, pen=penest, x.age=xx, lower=lower, upper=upper)
     
 }
 else out <- list(coefficients=x$coefficients, pen70=x$pen70.est, pen=penest, x.age=penout$x.age)
 
 if(print){
   cat("Coefficients: \n")
   print(x$coefficients)
   cat("\nPenetrance (%) by age 70: \n")
   print(x$pen70.est)
  }
 invisible(out)
}