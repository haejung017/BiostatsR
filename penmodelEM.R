penmodelEM <- function(parms, vbeta, data, design="pop", base.dist="Weibull", robust=FALSE, method="data", mode="dominant", q=0.02){

  agemin <- attr(data, "agemin")
  if(is.null(agemin)) stop("agemin is not found. Specify agemin for data by attr(data,\"agemin\") ")
  
  newdata <- carrierprobgeno(data, method=method, mode=mode, q=q)
  
    if(base.dist=="lognormal"){
  		if(parms[2] <= 0) stop("parms[2] has to be > 0")
		else
  	  	est1 <- nlm(loglik, c(parms[1], log(parms[2]), vbeta), data=data, design=design, base.dist=base.dist, agemin=agemin, hessian=TRUE)
  	} 
  	else if(any(parms<0)) stop("both parms have to be > 0")

  
  
  if(base.dist=="lognormal") theta = theta0 = c(parms[1],log(parms[2]), vbeta)
  else theta = theta0 = c(log(parms), vbeta)
  
  
  est0 <- est <- theta
  dd <- lval0 <- lval <- 1
  i <- 0
  lval <- loglikem(est, est0, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=FALSE)

  while(dd>0.00001){
    i <- i+1
    est0 <- est
    lval0 <- lval
    nlm.est <- nlm(loglikem, est0, theta0=est0, data=newdata, design=design, base.dist=base.dist, agemin=agemin, 
                   vec=FALSE, hessian=TRUE)
    lval <- nlm.est$minimum
    est <- nlm.est$estimate
    dd <- abs(lval0-lval)
    #dd <- abs(sum(est-est0))    
    #print(c(i, dd, lval, est))
  }
cat("Iterations = ", i, "\n")

  logLik <- -lval
  EST <- nlm.est$estimate
  Var <- try(solve(nlm.est$hessian), TRUE)

  if(!is.null(attr(Var,"class"))) stop("Model didn't converge.\n  Try again with different initial values")
  else{  
  	if(robust){
  	grad <- numericGradient(loglikem, nlm.est$estimate, theta0=nlm.est$estimate, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=TRUE)
  	Jscore <- t(grad)%*%grad
  	parms.cov <- Var%*%(Jscore)%*%Var
  	parms.se <- sqrt(diag(parms.cov))
  	}
  	else{   	parms.cov <- Var
  		  	parms.se <- sqrt(diag(parms.cov))
  	}
  }
  

  if(base.dist=="lognormal"){
  names(EST)<- names(parms.se)  <- c("lambda","log(rho)" , "beta.sex","beta.gene")
  rownames(parms.cov) <- colnames(parms.cov) <- c("lambda","log(rho)" , "beta.sex","beta.gene")
  }
  else{
  names(EST)<- names(parms.se)  <- c("log(lambda)","log(rho)" , "beta.sex","beta.gene")
  rownames(parms.cov) <- colnames(parms.cov) <- c("log(lambda)","log(rho)" , "beta.sex","beta.gene")
  }

  ageonset <- agemin:90
  
  p1 <- penf(nlm.est$estimate, ageonset, sex=1, mut=1, base.dist=base.dist, agemin=agemin)  
  p2 <- penf(nlm.est$estimate, ageonset, sex=0, mut=1, base.dist=base.dist, agemin=agemin)  
  p3 <- penf(nlm.est$estimate, ageonset, sex=1, mut=0, base.dist=base.dist, agemin=agemin)  
  p4 <- penf(nlm.est$estimate, ageonset, sex=0, mut=0, base.dist=base.dist, agemin=agemin)  
  
  pen.est <- pen.ci(nlm.est$estimate, parms.cov, age=70, base.dist=base.dist, agemin=agemin)
  pen70.est <- pen.est[1,]
  pen70.se <- pen.est[2,]
  pen70.ci <- pen.est[3:4,] 
  rownames(pen70.ci) <- c("lowerlimit", "upperlimit")
  
  
  out <- list(  coefficients=EST, varcov=parms.cov, se=parms.se,
 		pen70.est=pen70.est, pen70.se=pen70.se, pen70.ci=pen70.ci,ageonset=ageonset,  
        pen.maleCarr=p1, pen.femaleCarr=p2, pen.maleNoncarr=p3, pen.femaleNoncarr=p4,
        logLik=logLik)
  
  class(out) <- "penmodel"
  attr(out, "design") <- design
  attr(out, "base.dist") <- base.dist
  attr(out, "agemin") <- agemin
  attr(out, "data") <- data
  attr(out, "iterations") <- i
  attr(out, "robust") <- robust
  
invisible(out)

  }
