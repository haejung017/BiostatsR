# Carrier prob conditional on the phenotype
carrierprobpheno <- function(method="data", fit=NULL, data, mode="dominant", q=0.02)
{
  
  if(sum(is.na(data$mgene))==0) stop("Mutatioin carrier statuses are all known")
  if(method=="data"){
    carrp <- data$mgene
    
    cfam.id <- data$famID[data$proband==1 & data$mgene==1]
    nfam.id <- data$famID[data$proband==1 & data$mgene==0]
    i.cfam <- is.element(data$famID,cfam.id)
    i.nfam <- is.element(data$famID,nfam.id)
    
    for(g in unique(data$relation)){
      for(s in c(0,1)){
        for(d in c(0,1)){
          # carrier families
          carrp[i.cfam & is.na(data$mgene) & data$relation==g & data$gender==s & 
                  data$status==d] <- mean(data$mgene[i.cfam & !is.na(data$mgene) & data$relation==g & data$gender==s & data$status==d])
          # non-carrier famiiies
          carrp[i.nfam & is.na(data$mgene) & data$relation==g & data$gender==s & 
                  data$status==d] <- mean(data$mgene[i.nfam & !is.na(data$mgene) & data$relation==g & data$gender==s & data$status==d])
        }
      }
    }
  }  # close for method=="data"
  else if(method=="model"){
    if(is.null(fit)) stop("fit has to be specified.")
    theta <- fit$parms.est 
    base.dist <- attr(fit, "base.dist")
    agemin <- attr(fit, "agemin")

    beta.sex <- theta[3]
    beta.gen <- theta[4]
    xbeta1 <- beta.sex*data$gender+beta.gen*1
    xbeta0 <- beta.sex*data$gender+beta.gen*0
    y <- data$time-agemin
    delta <- data$status
    haz <-hazards(base.dist, y, theta)
    Haz <-cumhaz(base.dist, y, theta)
    haz1 <- haz*exp(xbeta1)
    haz0 <- haz*exp(xbeta0)
    Haz1 <- Haz*exp(xbeta1)
    Haz0 <- Haz*exp(xbeta0)
    p <- data$carrp.geno
    if(is.null(p)) {
      p <- carrierprobgeno(data=data, method="data", mode=mode, q=q)$carrp.geno
      data$carrp.geno <- p
    }
    p1 <- (haz1^delta)*exp(-Haz1)*p
    p0 <- (haz0^delta)*exp(-Haz0)*(1-p)
    carrp <- p1/(p1+p0)
    carrp[!is.na(data$mgene)] <- data$mgene[!is.na(data$mgene)]
  } # close for method=="model

  carrp[is.na(carrp)] <- 0
  data$carrp.pheno <- carrp
  
  return(data)

}