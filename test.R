for(n in c(500)){
  p <- 8; interest_lst <- 1:8; y_type <- "binary"; alpha <- -1; thres_w = 0.01; missing.rate <- 0.5
  variance <- F # If variance == T, the variance will be calculated

  if(variance == F){
    library(doParallel);
    closeAllConnections()
    n_cores <-  detectCores(all.tests = FALSE, logical = TRUE) # as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
    print(n_cores)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl, cores = n_cores)

    system.time(res <- foreach(index = 1:500, .packages = c("mgcv", "glmnet", "pROC", "MAVE"), .combine=rbind) %dopar%{
      set.seed(index)
      Simulation(p = p, n = n, interest_lst = interest_lst, thres_w = thres_w, y_type = y_type, intercept = T, variance = F, alpha = alpha, missing.rate = missing.rate)
    })
    stopCluster(cl)

  }else{ # bootstrap to calculate the variance

    library(doParallel);
    closeAllConnections()
    n_cores <-  detectCores(all.tests = FALSE, logical = TRUE) # as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
    print(n_cores)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl, cores = n_cores)

    system.time(foreach(index = 1:500, .packages = c("mgcv", "glmnet", "pROC", "MAVE")) %dopar%{
      set.seed(index)
      res.boot <- Simulation(p = p, n = n, interest_lst = interest_lst, thres_w = 0.01, y_type = y_type,intercept =T, variance = T, alpha = alpha, missing.rate = missing.rate)
      save(res.boot, file = paste0(n,"_", p,"_", y_type,"_",missing.rate,".RData"))
    })
    stopCluster(cl)
  }




}


