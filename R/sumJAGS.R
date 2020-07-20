sumJAGS <- function(rawJAGS, chunk.length = 50) {
  require(coda)
  require(stringr)
  require(testit)
  
  npars <- dim(rawJAGS$AA)[2]
  if(npars <= chunk.length) {
    Rhat <- gelman.diag(rawJAGS)$psrf[, 2]
    neff <- effectiveSize(rawJAGS)
  } else { # Avoids "Error in chol.default(W) :" with large models
    Rhat <- neff <- c()
    nchunks <- ceiling(npars / chunk.length)
    for(chnk in 1:nchunks) {
      st <- (chnk*chunk.length) - chunk.length + 1
      end <- ifelse(chnk == nchunks, npars, chnk*chunk.length)
      if(has_error(gelman.diag(rawJAGS[, st:end]), silent = T)) {
        for(i in st:end) Rhat <- c(Rhat, gelman.diag(rawJAGS[, i])$psrf[, 2])
      } else {
        Rhat <- c(Rhat, gelman.diag(rawJAGS[, st:end])$psrf[, 2])
      }
      neff <- c(neff, effectiveSize(rawJAGS[, st:end]))
    }
  }
  
  sims.list <- simsList(rawJAGS)
  cols <- c("mean", "SD", "q025", "q250", "q500", "q750", "q975", "Rhat", "neff")
  rows <- names(Rhat)
  sumTab <- matrix(NA, nrow = length(rows), ncol = length(cols),
                   dimnames = list(rows, cols))
  sumrow.fn <- function(x) c(mean(x), sd(x),
                             quantile(x, prob = 0.025, type = 8),
                             quantile(x, prob = 0.25, type = 8),
                             quantile(x, prob = 0.5, type = 8),
                             quantile(x, prob = 0.75, type = 8),
                             quantile(x, prob = 0.975, type = 8))

  for(i in 1:nrow(sumTab))
    if(!str_detect(rows[i], "\\["))
      sumTab[rows[i], cols[1:(length(cols) - 2)]] <- sumrow.fn(sims.list[[rows[i]]])
  ind.p.arrays <- which(str_detect(rows, "\\["))
  if(length(ind.p.arrays) > 0) {
    p.arrays <- unique(str_split(rows[ind.p.arrays], "\\[", simplify = T)[,1])
    for(p in p.arrays) {
      p.ests <- sims.list[[p]]
      p.mean <- apply(p.ests, 2:length(dim(p.ests)), mean)
      p.sd <- apply(p.ests, 2:length(dim(p.ests)), sd)
      p.025 <- apply(p.ests, 2:length(dim(p.ests)), function(x) quantile(x, prob = 0.025, type = 8))
      p.250 <- apply(p.ests, 2:length(dim(p.ests)), function(x) quantile(x, prob = 0.25, type = 8))
      p.500 <- apply(p.ests, 2:length(dim(p.ests)), function(x) quantile(x, prob = 0.5, type = 8))
      p.750 <- apply(p.ests, 2:length(dim(p.ests)), function(x) quantile(x, prob = 0.75, type = 8))
      p.975 <- apply(p.ests, 2:length(dim(p.ests)), function(x) quantile(x, prob = 0.975, type = 8))
      p.est.dims.chr <- rows[which(str_sub(rows, 1, nchar(p)) == p)]
      if(any(!str_detect(p.est.dims.chr, "\\[")))
        p.est.dims.chr <- p.est.dims.chr[-which(!str_detect(p.est.dims.chr, "\\["))]
      if(any(str_sub(p.est.dims.chr, 1, nchar(p) + 1) != str_c(p, "[")))
        p.est.dims.chr <- p.est.dims.chr[-which(str_sub(p.est.dims.chr, 1, nchar(p) + 1) != str_c(p, "["))]
      p.est.dims.int <- str_split(p.est.dims.chr, "\\[", simplify = T)[,2]
      p.est.dims.int <- str_split(p.est.dims.int, "\\]", simplify = T)[,1]
      if(any(str_detect(p.est.dims.int, ","))) {
        p.est.dims <- str_split(p.est.dims.int, ",", simplify = T)
        dim.p.est.dims <- dim(p.est.dims)
        p.est.dims <- array(as.integer(p.est.dims), dim = dim.p.est.dims)
      } else {
        p.est.dims <- as.integer(p.est.dims.int)
      }
      sumTab[p.est.dims.chr, cols[1:(length(cols) - 2)]] <- cbind(
        as.numeric(p.mean), as.numeric(p.sd), as.numeric(p.025),
        as.numeric(p.250), as.numeric(p.500), as.numeric(p.750),
        as.numeric(p.975))
    }
  }
  sumTab[, "Rhat"] <- Rhat
  sumTab[, "neff"] <- neff
  
  out <- list(sims.array = rawJAGS, sims.list = sims.list, summary = sumTab)
}
