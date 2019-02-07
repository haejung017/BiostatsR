print.summary.penmodel <- function (x, digits = max(3, getOption("digits") - 3), 
                                    signif.stars=TRUE, ...) 
{
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Coefficients: \n")

  
  if(signif.stars){
    Signif <- symnum(x$coefficients[,4], corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " "))
    
    coefstars <- cbind(data.frame(x$coefficients), format(Signif) )
    colnames(coefstars) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
    print(coefstars)
    cat("Signif. codes:   0 '***'  0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
    x$coefficients <- coefstars 
  }
  else print(x$coefficients)
    
  cat("\nPenetrance (%) by age 70: \n")
  print(x$pen70[1:2,])
  cat("\n95% Confidence intervals on the penetrances: \n")
  print(x$pen70[3:4,])
  
  invisible(x)
  
}
