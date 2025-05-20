
#Function to filter out outliers log-likelihoods using IQR
filt <- function(x){
  x1 <- x[!is.infinite(x) & !is.nan(x)]
  wid <- 1.5*IQR(x1,na.rm=TRUE)
  low <- quantile(x1,0.25,na.rm=TRUE) - wid
  high <- quantile(x1,0.75,na.rm=TRUE) + wid
  for(j in 1:length(x)){
    xj <- x[j]
    if(xj > high || xj < low || is.nan(xj) == TRUE || is.infinite(xj) == TRUE){
      x[j] = NA
    }
  }
  return(x)
}