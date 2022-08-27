## Sample generation ####
sample.generation <- function(p = 8, n = 200, y_type = c("continuous", "binary"), alpha = -1, missing.rate = 0.5){ # missing.rate = "medium"

  if(!(y_type %in% c("continuous", "binary")) ) stop("y_type should be either continuous or binary")
  if(p <= 4) stop("p should be larger than 4")
  ## R
  R <- rbinom(n, 1, 1-missing.rate)

  ## Covariates
  X <- matrix(NA, nrow = n, ncol = p)
  var_p <- diag(p)
  # var_p2 <- 1.5*0.9^(0:7)
  # var_p2 <- matrix(rbind(var_p2, do.call("rbind", lapply(1:7, function(idx){
  #   1.5 * 0.9^c(idx:1, 0:(7-idx))
  # }))), ncol=8)
  # var_p2 <- diag(rep(c(2,1), each=4))
  var_p2 <- var_p*1.5

  X[R==0,] <- mgcv::rmvn(sum(R==0), rep(0,p), V=var_p); # D1

  gamma <- .3
  mixture <- rbinom(sum(R==1), 1, gamma)
  X[R==1,][mixture==1,] <- mgcv::rmvn(sum(mixture==1), rep(0,p), V=var_p) # D1
  X[R==1,][mixture==0,] <- mgcv::rmvn(sum(mixture==0), rep(1,p), V=var_p2); # D2 var_p*1.5; diag(rep(c(2,1), each=4))


  colnames(X) <- paste0("X", 1:p)

  ## True beta
  if(y_type == 'continuous'){
    beta0 <- c(1,1,-1,-1, rep(0, times=p-4))
  }else if(y_type == 'binary'){
    beta0 <- 0.5*c(1,1,-1,-1, rep(0, times=p-4))
  }

  ## Outcomes Z and Y2
  Z <- X %*% beta0/2 + rowSums(abs(X[,5:p]))/(p-4) + rnorm(n, sd = 1)

  # Y2
  Y2 <- (1-alpha) * Z + X %*% beta0/2 + rnorm(n, sd = 1)


  ## function h
  if(y_type == "continuous"){
    y <- Y2-Z
    impute_true <- -alpha*Z + X%*% beta0
  }else if(y_type == "binary"){
    y <- (Y2 >= Z)*1
    impute_true <- pnorm(alpha*Z - X%*%beta0, mean = 0, sd = 1, lower.tail = F)

  }
  y[R==0] <- NA

  list(p = p, n = n, missing.rate=mean(R), X = X, Z = Z, R = R, y=y, impute_true = impute_true, missing_true = rep(1-missing.rate, n)) # Y2 = Y2,
}

### Simulation ####

##### SIMULATION SETTING #####

# load estimate, debias, sample.generation

# library(mgcv); library(glmnet); library(parallel); library(MAVE)


### Looking at the beta ###
Simulation <- function(p = 8, n = 200, y_type = 'continuous', interest_lst = 1:4, intercept = F, thres_w = 0.01, variance = F, alpha=0, missing.rate = 0.5){
  samples <- sample.generation(p = p, n= n, y_type = y_type, alpha = alpha, missing.rate =missing.rate)

  X <- samples$X; Z <- samples$Z; y <- samples$y # R <- samples$R;
  nuisance_true <- list(impute = samples$impute_true, missingness = samples$missing_true)


  results_boot <- NULL; results_beta <- NULL;
  if(variance == F){
    #     results_true <- missInfer(samples = samples_re, nuisance_true = nuisance_true, interest_lst = interest_lst,
    #                                  intercept = intercept, thres_w =  thres_w, variance = F,  conditioning = "", pi.method = "true", dim.reduction = "")
    #     results_true <- unlist(results_true)
    #     names(results_true) <- paste0(names(results_true), "_", rep(interest_lst, each = 2))


    results_glmnet <- missInfer(X = X, y = y, interest_var = 1:p,
                                   intercept = intercept, thres_w =  thres_w, pi.method = "glmnet", dim.reduction = "onlyQ")
    names(results_glmnet) <- paste0("glmnet_",names(results_glmnet))

    results_beta_x_sep <- missInfer(X = X, y = y, interest_var = 1:p,
                                    intercept = intercept, thres_w =  thres_w, pi.method = "kernel", dim.reduction = "separate")
    names(results_beta_x_sep ) <- paste0("x_sep_",names(results_beta_x_sep))

    results_beta_x <- missInfer(X = X, y = y, interest_var = 1:p,
                                intercept = intercept, thres_w =  thres_w, pi.method = "kernel", dim.reduction = "onlyQ")
    names(results_beta_x ) <- paste0("x_",names(results_beta_x))

    results_beta_xy <- missInfer(X = X, Z = Z, y = y, interest_var = 1:p,
                                 intercept = intercept, thres_w =  thres_w, pi.method = "kernel", dim.reduction = "onlyQ")
    names(results_beta_xy ) <- paste0("xy_",names(results_beta_xy))



    results_beta <- c(
      #                      results_true, true_prop.pihat = sum(samples_re$missing_true < 0.005)/n,
      # results_beta_x_sep, results_glmnet, glmnet_prop.pihat = sum(predict(mod_missingness, X, type='response', s = "lambda.min") <0.005)/n,
      results_beta_x_sep, results_glmnet,
      results_beta_x, results_beta_xy, missing.rate= 1-mean(samples$R), frac = mean(samples$y, na.rm = T) #,
      # supp.pi.X = 8, # ifelse(cv.mave.pi.X$dim.min > 0, cv.mave.pi.X$dim.min , p),
      # supp.pi.glmnet = sum(supp.pi.glmnet!=0),
      # supp.Q.X = 5, # ifelse(cv.mave.X$dim.min > 0, cv.mave.X$dim.min , p),
      # supp.Q.ZX = ifelse(cv.mave.ZX$dim.min > 0, cv.mave.ZX$dim.min , p)
      )

    # results_beta <- results_true
  }else if(variance == T){
    no_boot <- 500

    # results_boot <- foreach(boot = 1:no_boot, .packages = c("mgcv", "glmnet", "MAVE"), .combine = rbind) %dopar%{
    results_boot <- do.call("rbind", lapply(1:no_boot, function(boot){
      print(paste0("boot number: ",boot))
      set.seed(boot)

      err_val <- 0 ; count <- 0
      while(err_val == 0 ){
        boot_idx <- sample(n, n, replace = T)

        samples_boot <- list(n = n, p = p, X = samples_re$X[boot_idx,], Z = samples_re$Z[boot_idx],
                             Y2 = samples_re$Y2[boot_idx],  R = samples_re$R[boot_idx], y=samples_re$y[boot_idx],
                             impute_true = samples_re$impute_true[boot_idx], missingness_true = samples_re$missing_true[boot_idx])
        nuisance_true_boot <- list(impute = samples_re$impute_true[boot_idx], missingness = samples_re$missing_true[boot_idx])

        tryCatch({
          results_glmnet_boot <- missInfer(samples = samples_boot,interest_lst = interest_lst,
                                              intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "X", pi.method = "glmnet", dim.reduction = "onlyQ",
                                              supp.Q = supp.Q.X, supp.pi = supp.pi.glmnet)

          results_beta_x_sep_boot <- missInfer(samples = samples_boot, interest_lst = interest_lst,
                                                  intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "X", pi.method = 'kernel', dim.reduction = "separate",
                                                  supp.Q = supp.Q.X, supp.pi = supp.pi.X)

          results_beta_x_boot <- missInfer(samples = samples_boot, interest_lst = interest_lst,
                                              intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "X", pi.method = 'kernel', dim.reduction = "onlyQ",
                                              supp.Q = supp.Q.X, supp.pi = supp.Q.X)

          results_beta_xy_boot <- missInfer(samples = samples_boot,  interest_lst = interest_lst,
                                               intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "ZX", pi.method = 'kernel', dim.reduction = 'onlyQ',
                                               supp.Q = supp.Q.ZX, supp.pi = supp.Q.ZX)

          err_val <- 1
        }, error = function(err){
          err_val <- 0
          if(count >= 3){
            message(" Here is the original error message:")
            message(err)
            return(NA)
          }
          count <- count + 1
        })
      }

      results_glmnet_boot <- unlist(results_glmnet_boot)
      names(results_glmnet_boot) <- paste0("glmnet_",names(results_glmnet_boot), "_", rep(interest_lst, each = 2))

      results_beta_x_sep_boot <- unlist(results_beta_x_sep_boot)
      names(results_beta_x_sep_boot) <- paste0("x_sep_",names(results_beta_x_sep_boot), "_", rep(interest_lst, each=2+length(thres_w)))

      results_beta_x_boot <- unlist(results_beta_x_boot)
      names(results_beta_x_boot) <- paste0("x_",names(results_beta_x_boot), "_", rep(interest_lst, each=2+length(thres_w)))

      results_beta_xy_boot <- unlist(results_beta_xy_boot)
      names(results_beta_xy_boot) <- paste0("xy_",names(results_beta_xy_boot), "_", rep(interest_lst, each=2+length(thres_w)))

      c(results_beta_x_sep_boot, results_glmnet_boot, results_beta_x_boot, results_beta_xy_boot)


    }))

  }



  return(rbind(results_beta, results_boot) )

}





# for(n in c(500)){ # ,350, 500 # 500, 650
#   sample.size <- n; p <- 8; interest_lst <- 1:8; y_type <- "binary"; variance <- F; alpha <- -1
#   missing.rate <- 0.5
#   # sampling_lower <- AAA+1; sampling_upper <- BBB
#
#
#
#   if(variance == F){
#     library(doParallel);
#     closeAllConnections()
#     n_cores <-  detectCores(all.tests = FALSE, logical = TRUE) # as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
#     print(n_cores)
#     cl <- makeCluster(n_cores)
#     registerDoParallel(cl, cores = n_cores)
#
#     system.time(res <- foreach(index = 1:1, .packages = c("mgcv", "glmnet", "pROC", "MAVE"), .combine=rbind) %dopar%{
#
#
#       set.seed(index)
#       Simulation(p = p, n = n, interest_lst = interest_lst, thres_w = .01, y_type = y_type, intercept = T, variance = F, alpha = alpha, missing.rate = missing.rate)
#
#       # simulation.result
#
#     })
#     # save(res, file=paste0("Simulations/Inference/boot_120821/results/Sim_",sample.size,"_",p,"_", y_type,"_",alpha,"_",missing.rate,"_120821.RData"))
#     stopCluster(cl)
#   }else{
#
#     lapply(BBB:(BBB), function(AAA){
#       set.seed(AAA)
#       res.boot <- Simulation(p = p, n = n, interest_lst = interest_lst, thres_w = 0.01, y_type = y_type,intercept =T, variance = T, alpha = alpha, missing.rate = missing.rate)
#       # res.boot
#       # })
#       save(res.boot, file=paste0("Simulations/Inference/boot_120821/results/n",n,"_p",p,"_", y_type,"_90/Sim_",sample.size,"_",p,"_", y_type,"_",alpha, "_boot_", AAA,"_",missing.rate,"_120821.RData"))
#     })
#
#
#   }
#   #
#
#
#
# }
#
#
