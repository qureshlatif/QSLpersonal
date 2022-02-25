Impute_missing_covs_rf <- function(dat, v.fill, v.inform) {
  require(randomForest)
  require(dplyr)
  require(stringr)

  # Define data objects and convert to data frames as needed.
  dat.mat <- is.matrix(dat)
  if(dat.mat) dat <- data.frame(dat)

  ind.missing <- which(is.na(dat[,v.fill]))
  ind.known <- which(!is.na(dat[,v.fill]))
  rf <- randomForest(as.formula(str_c(v.fill, "~",
                                      str_c(v.inform[which(v.inform != v.fill)],
                                            collapse = "+"))),
                     data = (dat %>% slice(ind.known)))
  dat[ind.missing, v.fill] <- predict(rf, newdata = (dat %>% slice(ind.missing)))
  if(dat.mat) dat <- data.matrix(dat)
  return(dat)
}
