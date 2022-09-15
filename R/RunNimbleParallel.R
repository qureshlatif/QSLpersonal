RunNimbleParallel <-
  function(model, inits, data, constants, parameters,
           nc = 2, ni = 2000, nb = 0.5, nt = 10, mod.nam = "mod",
           max.samples.saved = 10000, rtrn.model = F,
           Rht.required = 1.1, neff.required = 100) {
    #~~~~~Parallel processing code (probably won't work in Windows)~~~~~#
    require(nimble)
    require(parallel)
    require(coda)
    require(mcmcOutput)
    set.seed(1)
    cl<-makeCluster(nc, timeout = 5184000)
    clusterExport(cl, c("model", "inits", "data", "constants", "parameters",
                        "ni", "nt"))
    for (j in seq_along(cl)) {
      set.seed(j)
      init <<- inits()
      clusterExport(cl[j], "init")
    }
    out1 <- clusterEvalQ(cl, {
      library(nimble)
      library(coda)
      model <- nimbleModel(code = model, name = "model",
                           constants = constants, data = data,
                           inits = init)
      Cmodel <- compileNimble(model)
      modelConf <- configureMCMC(model, thin = nt)
      modelConf$addMonitors(parameters)
      modelMCMC <- buildMCMC(modelConf)
      CmodelMCMC <- compileNimble(modelMCMC, project = model)
      CmodelMCMC$run(ni, reset = FALSE)
      return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
    })
    
    ## Check convergence ##
    out2 <- out1
    ni.saved <- nrow(out2[[1]])
    for(chn in 1:nc) { # nc must be > 1
      out2[[chn]] <- out2[[chn]][(round(ni.saved * nb)+1):ni.saved,]
    }
    out.mcmc <- coda::as.mcmc.list(lapply(out2, coda::as.mcmc))
    
    mod <- mcmcOutput(out.mcmc)
    sumTab <- summary(mod, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = F)
    sumTab <- sumTab %>%
      as_tibble() %>%
      mutate(Parameter = row.names(sumTab)) %>%
      select(Parameter, mean:f)
    
    ind.Rht <- which(!str_detect(sumTab$Parameter, "test.n") &
                       !str_detect(sumTab$Parameter, "M.save") &
                       !str_detect(sumTab$Parameter, "pXtheta") &
                       !str_detect(sumTab$Parameter, "Ind") &
                       !str_detect(sumTab$Parameter, "ind"))
    mxRht <- sumTab %>% slice(ind.Rht) %>% pull(Rhat) %>% max(na.rm = T)
    mn.neff <- sumTab %>% slice(ind.Rht) %>% pull(n.eff) %>% min(na.rm = T)

    mod <- list(mcmcOutput = mod, summary = sumTab)
    R.utils::saveObject(mod, mod.nam) # If running all in one.
    
    ## If has not converged, continue sampling
    if(round(mxRht, digits = 1) > Rht.required | mn.neff < neff.required) {
      n.runs <- 1
      R.utils::saveObject(out1, str_c(mod.nam, "_chunk", n.runs)) # Save samples from previous run to drive.
    }
    while(round(mxRht, digits = 1) > Rht.required | mn.neff < neff.required) {
      n.runs <- n.runs + 1
      print(str_c("Run = ", n.runs))
      
      out2 <- clusterEvalQ(cl, {
        CmodelMCMC$run(ni, reset = FALSE, resetMV = TRUE) # Resume sampling.
        return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
        gc(verbose = F)
      })
      R.utils::saveObject(out2, str_c(mod.nam, "_chunk", n.runs)) # Save samples from previous run to drive.
      
      ni2 <- round(((ni / nt) * n.runs * nc) * nb) # Anticipated number of samples to save (assuming half discarded as burn-in).
      if(ni2 > max.samples.saved) {
        nt2 <- round(1 / (max.samples.saved / ni2)) # Set additional thinning so that saved iterations don't exceed (by too much) max.samples.saved (specified by user).
      } else {
        nt2 <- 1
      }
      
      # Reassemble chain from chunks and apply additional thinning.
      out1 <- R.utils::loadObject(str_c(mod.nam, "_chunk", 1))
      ni.saved <- nrow(out1[[1]])
      for(chn in 1:nc) { # nc must be > 1
        out1[[chn]] <- out1[[chn]][seq(2, ni.saved, by = nt2),] # Starting at 2 because first iteration is NA for some reason.
      }
      for(r in 2:n.runs) {
        out.r <- R.utils::loadObject(str_c(mod.nam, "_chunk", r))
        ni.saved <- nrow(out.r[[1]])
        for(chn in 1:nc) {
          out.r[[chn]] <- out.r[[chn]][seq(1, ni.saved, by = nt2),]
          out1[[chn]] <- rbind(out1[[chn]], out.r[[chn]])
        }
      }
      
      # Discard first half of samples as burn-in
      out3 <- out1
      ni.saved <- nrow(out3[[1]])
      for(chn in 1:nc) {
        out3[[chn]] <- out3[[chn]][(round(ni.saved/2)+1):ni.saved,]
      }
      out.mcmc.update <- coda::as.mcmc.list(lapply(out3, coda::as.mcmc))
      
      mod <- mcmcOutput(out.mcmc.update)
      sumTab <- summary(mod, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = F)
      sumTab <- sumTab %>%
        as_tibble() %>%
        mutate(Parameter = row.names(sumTab)) %>%
        select(Parameter, mean:f)
      mxRht <- sumTab %>% slice(ind.Rht) %>% pull(Rhat) %>% max(na.rm = T)
      gc(verbose = F)
      
      mod <- list(mcmcOutput = mod, summary = sumTab)
      R.utils::saveObject(mod, mod.nam) # If running all in one.
    }
    for(r in 1:n.runs) file.remove(str_c(mod.nam, "_chunk", r))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if(rtrn.model) return(mod)
  }
  