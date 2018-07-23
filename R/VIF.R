VIF <- function(data, vars = names(data), returnModels = F) {
  vif <- numeric(length = length(vars))
  names(vif) <- vars
  dat <- data[vars]
  mods <- vector("list", length = length(vars))
  names(mods) <- vars
  for(i in 1:length(vars)) {
    mod <- lm(as.formula(paste0(vars[i], "~", paste(vars[-i], collapse = "+"))), data = dat)
    vif[i] <- 1 / (1 - summary(mod)$r.squared)
    mods[[i]] <- mod
  }
  ifelse(returnModels == T,
         return(list(vif = vif, mods = mods)),
         return(vif))
}
