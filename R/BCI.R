BCI <- function(x, ndig = 2, BCIpercent = 95, flag.sig = T) {
  plo <- 0 + (1 - BCIpercent/100)/2
  phi <- 1 - (1 - BCIpercent/100)/2
  md <- median(x)
  lo <- quantile(x, prob = plo, type = 8)
  hi <- quantile(x, prob = phi, type = 8)
  x.sum <- ifelse((lo > 0 | hi < 0) & flag.sig,
                  str_c(round(md, digits = ndig),
                        " (",
                        round(lo, digits = ndig),
                        ",",
                        round(hi, digits = ndig),
                        ")*"),
                  str_c(round(md, digits = ndig),
                        " (",
                        round(lo, digits = ndig),
                        ",",
                        round(hi, digits = ndig),
                        ")"))
  return(x.sum)
}                          
