



fun_debias <- function(samples, nuisance_true = NULL, interest_var = 1, intercept = F, thres_w = 0, variance=F, conditioning = c("X", "ZX"),
                         pi.method = c("true", "kernel", "glmnet"), dim.reduction = c("separate", "onlyQ"),
                         supp.pi = NULL, supp.Q = NULL, beta_initial = NULL){

  #### variable assignments #####
  X <- samples$X; Z <- samples$Z; R <- samples$R; y <- samples$y
  n <- nrow(X); p <- ncol(X);
  if(length(unique(y)[!is.na(unique(y))])==2 ){
    y_type <- "binary"
  }else { y_type <- "continuous"}

  if(intercept){
    interest_var <- c(0, interest_var)
  }

  #### model for missingness ####
  if(pi.method == "glmnet"){
    if(intercept){
      missingness_glmnet <- b_1st(cbind(1, X), beta_tilde = supp.pi, y_type = "binary")
    }else{
      supp.pi <- supp.pi[-1]
      missingness_glmnet <- b_1st(X, beta_tilde = supp.pi, y_type = "binary")
    }

    missingness_glmnet <- truncated(missingness_glmnet, thres = thres_w)
  }


  if(conditioning == "X"){
    combined_dataset <- X
  }else if(conditioning == "ZX"){
    combined_dataset <- cbind(Z=Z, X)
  }


  if(pi.method %in% c("glmnet", "kernel")){
    ### STEP 1 ####
    # dim.reduction.Q <- mave(y ~ ., data = data.frame(y=y[R==1], combined_dataset[R==1,]), method =
    #                           "KSIR") #"KSIR"
    # cv.mave <- mave.dim(dim.reduction.Q)


    ### STEP 2 ####
    # print(combined_dataset)


    if(is.null(beta_initial)){
      if(!is.null(supp.Q)){
        imputed_entire <- ks(xx=as.matrix(combined_dataset[R==1,]) %*% supp.Q, yy = y[R==1], xx.test = as.matrix(combined_dataset) %*% supp.Q)
        #### here
      }else{
        imputed_entire <- ks(xx=as.matrix(combined_dataset[R==1,]), yy = y[R==1], xx.test = as.matrix(combined_dataset))
      }

      beta_initial <- get_initial(samples = samples, impute = imputed_entire, intercept = intercept)
      # beta_initial <- coef(mod_impute)
    }


  }


  sim_results <- lapply(interest_var, function(interest_idx){



    ### STEP 3 ####

    if(intercept){ ### Still working

      if(interest_idx == 0){
        mod_w <- glmnet::glmnet(x =  X, y = rep(1, nrow(X)), intercept=FALSE, weights = b_2nd(cbind(1, X), beta_initial, y_type), lambda = 0)
        w_est <- coef(mod_w)[-1,]
      }else{
        mod_w <- glmnet::glmnet(x =  X[, -interest_idx], y = X[,interest_idx], intercept=intercept, weights = b_2nd(cbind(1, X), beta_initial, y_type), lambda = 0)
        w_est <- coef(mod_w)
      }

      if(interest_idx ==0){
        w_tilde <- c(1, -w_est[1:p])
      } else if(interest_idx >=1 & interest_idx <p){
        w_tilde <- c(-w_est[1:(interest_idx)], 1, -w_est[(interest_idx):(p-1)])
      } else{
        w_tilde <- c(-w_est[1:p], 1)
      }


    }else{


      mod_w <- glmnet::glmnet(x =  X[, -interest_idx], y = X[,interest_idx], intercept=intercept, weights = b_2nd(X, beta_initial[-1], y_type), lambda = 0) #beta_initial[-1]
      w_est <- coef(mod_w)[-1,]

      if(interest_idx ==1){
        w_tilde <- c(1, -w_est)
      } else if(interest_idx >=2 & interest_idx <p){
        w_tilde <- c(-w_est[1:(interest_idx-1)], 1, -w_est[(interest_idx):(p-1)])
      } else{
        w_tilde <- c(-w_est, 1)
      }
    }

    # }

    # set.seed(1823)
    if(pi.method %in% c("glmnet", "kernel")){
      ### STEP 4 ####
      sampleSplitIndex <- rep(F,n)
      sampleT <- sample(which(R == 1), size = round(sum(R==1)/2))
      sampleF <- sample(which(R == 0), size = round(sum(R==0)/2))
      sampleSplitIndex[c(sampleT, sampleF)] <- T

      ### STEP 5 ####

      nuisance <- lapply(c(T,F), function(sample.div){

        X_sub <- X[sampleSplitIndex==sample.div, ]
        combined_sub <- combined_dataset[sampleSplitIndex==sample.div,];
        y_sub <- y[sampleSplitIndex==sample.div]; R_sub <- R[sampleSplitIndex==sample.div];

        X_opp <- X[sampleSplitIndex!=sample.div, ]
        combined_opp <- combined_dataset[sampleSplitIndex!=sample.div,];
        y_opp <- y[sampleSplitIndex!=sample.div]; R_opp <- R[sampleSplitIndex!=sample.div];

        if(NCOL(Z)==1){
          Z_sub <- Z[sampleSplitIndex==sample.div]
          Z_opp <- Z[sampleSplitIndex!=sample.div]
        }else{
          Z_sub <- Z[sampleSplitIndex==sample.div,]
          Z_opp <- Z[sampleSplitIndex!=sample.div,]
        }


        # print(combined_sub[R_sub==1,])

        if(!is.null(supp.Q)){
          impute_sub <- ks(xx=as.matrix(combined_sub[R_sub==1,])%*% supp.Q, yy = y_sub[R_sub==1], xx.test = as.matrix(combined_opp) %*% supp.Q)
        }else{
          impute_sub <- ks(xx=as.matrix(combined_sub[R_sub==1,]), yy = y_sub[R_sub==1], xx.test = as.matrix(combined_opp))
        }

        if(pi.method == "kernel"){


          if(dim.reduction=="separate"){
            # dim.reduction.pi <- mave(R ~ ., data = data.frame(R=R_sub, combined_sub), method ="KSIR") #KSIR
            # cv.mave.pi <- mave.dim(dim.reduction.pi)
            if(!is.null(supp.pi)){
              # supp.pi <- as.matrix(dim.reduction.pi$dir[[cv.mave.pi$dim.min]])
              # supp.pi <- as.matrix(dim.reduction.pi$dir[[length(dim.reduction.pi$dir)]])
              pi_conditional <- ks(xx = as.matrix(combined_sub) %*% supp.pi, yy=R_sub, xx.test = as.matrix(combined_opp) %*% supp.pi)
            }else{
              # supp.pi <- combined_sub
              pi_conditional <- ks(xx = as.matrix(combined_sub), yy=R_sub, xx.test = as.matrix(combined_opp))
            }

            pi_truncated <- lapply(thres_w, function(thres){
              tmp <- pi_conditional
              thres_below <- (tmp < thres)
              tmp[thres_below  ] <- thres
              #               list(tmp = tmp, extremely_small = sum(thres_below))
              tmp
            })

          }else{
            supp.pi <- supp.Q

            if(!is.null(supp.pi)){
              if(intercept){
                numerator <- ks.missing(xx = as.matrix(combined_sub[R_sub==1, ]) %*% supp.pi, yy=cbind(1, X_sub[R_sub == 1,]) %*% w_tilde, samp_size = sum(R_sub==1), xx.test = as.matrix(combined_opp) %*% supp.pi)
                denominator <- ks.missing(xx = as.matrix(combined_sub[R_sub==0,]) %*% supp.pi, yy=cbind(1, X_sub[R_sub == 0,]) %*% w_tilde,samp_size =  sum(R_sub==0), xx.test = as.matrix(combined_opp) %*% supp.pi)
              }else{
                numerator <- ks.missing(xx = as.matrix(combined_sub[R_sub==1, ]) %*% supp.pi, yy=X_sub[R_sub == 1,] %*% w_tilde, samp_size = sum(R_sub==1), xx.test = as.matrix(combined_opp) %*% supp.pi)
                denominator <- ks.missing(xx = as.matrix(combined_sub[R_sub==0,]) %*% supp.pi, yy=X_sub[R_sub == 0,] %*% w_tilde,samp_size =  sum(R_sub==0), xx.test = as.matrix(combined_opp) %*% supp.pi)
              }

            }else{
              if(intercept){
                numerator <- ks.missing(xx = as.matrix(combined_sub[R_sub==1, ]) , yy=cbind(1, X_sub[R_sub == 1,]) %*% w_tilde, samp_size = sum(R_sub==1), xx.test = as.matrix(combined_opp))
                denominator <- ks.missing(xx = as.matrix(combined_sub[R_sub==0,]), yy=cbind(1, X_sub[R_sub == 0,]) %*% w_tilde, samp_size = sum(R_sub==0), xx.test = as.matrix(combined_opp))
              }else{
                numerator <- ks.missing(xx = as.matrix(combined_sub[R_sub==1, ]) , yy=X_sub[R_sub == 1,] %*% w_tilde, samp_size = sum(R_sub==1), xx.test = as.matrix(combined_opp))
                denominator <- ks.missing(xx = as.matrix(combined_sub[R_sub==0,]), yy=X_sub[R_sub == 0,] %*% w_tilde, samp_size = sum(R_sub==0), xx.test = as.matrix(combined_opp))
              }
            }



            w_conditional <- numerator/denominator

            rho <- mean(R_opp)


            pi_conditional <- w_conditional * rho/(rho * w_conditional + 1 - rho)

            pi_truncated <- lapply(thres_w, function(thres){
              tmp <- pi_conditional
              tmp[(tmp > 0) & (tmp < thres) ] <- thres
              tmp[(tmp < 0) & (tmp > -thres) ] <- -thres
              tmp
            })


          }



          list(lapply(1:length(thres_w), function(idx){

            list(impute = impute_sub, missingness = pi_truncated[[idx]],
                 extremely_small_pi = ifelse(variance==F, sum(abs(pi_truncated[[idx]]) < thres_w[idx]), NA) )

          }))

        }else{

          list(impute = impute_sub, missingness =  missingness_glmnet[sampleSplitIndex!=sample.div,])
        }





      })
    }

    # print(c(S_bar, I_bar));

    if(pi.method %in% c("glmnet", "kernel")){
      I_bar <- fun_I(X=X, interest_idx = interest_idx, beta_tilde = beta_initial, w_hat = w_est, intercept=intercept, y_type = y_type)
      S_bar <- lapply(c(T,F), function(sample.div){
        samples_sub <- make_sub(samples, sampleSplitIndex, sample.div)

        # S_ks <-
        if(pi.method =='kernel'){

          lapply(1:length(thres_w), function(idx){
            nuisance_sample <- nuisance[[1+sample.div]][[1]][[idx]]
            loss <- loss_1st_boot(samples_sub, nuisance_sample, beta_est = beta_initial, y_type = y_type, intercept = intercept)
            S <- t(loss) %*% w_tilde / sum(sampleSplitIndex == sample.div)
            S
          })
        }else{
          nuisance_sample <- nuisance[[1+sample.div]]
          loss <- loss_1st_boot(samples_sub, nuisance_sample, beta_est = beta_initial, y_type = y_type, intercept = intercept)
          S <- as.numeric(t(loss) %*% w_tilde / sum(sampleSplitIndex == sample.div))
          S
        }

      })
    }
    # print(c(S_bar, I_bar));

    if(pi.method =="true" & variance == F){
      nuisance_true$missingness[nuisance_true$missingness < thres_w] <- thres_w
      I_bar_true <- fun_I(X=X, interest_idx = interest_idx, beta_tilde = beta_initial, w_hat = w_est, intercept = F, y_type = y_type)
      S_bar_true <- t(loss_1st_boot(samples, nuisance_true, beta_est =  beta_initial, y_type = y_type, intercept = F)) %*% w_tilde / n

    }else{I_bar_true <- NULL; S_bar_true <- NULL}


    if(pi.method == "kernel"){
      if(variance ==F){
        prop.pihat <- sapply(1:length(thres_w), function(idx){ sum(sapply(lapply(lapply(nuisance, '[[',1), '[[',idx),'[[',3))/n})
      }else{prop.pihat <- NA}
      c(
        # initial = beta_initial[interest_idx+1,],
        debiased = beta_initial[interest_idx+1,]-
          sapply(1:length(thres_w), function(idx){colMeans(do.call('rbind', lapply(S_bar, "[[", idx)))}) / I_bar)
    }else if(pi.method == "glmnet"){
      c(
        # initial = beta_initial[interest_idx+1,],
        debiased = beta_initial[interest_idx+1,] -
          mean(unlist(S_bar)) / I_bar)
    }else if(pi.method == "true"){
      c(
        # true_initial = beta_initial[interest_idx+1,],
        true_debiased =  beta_initial[interest_idx+1,] - S_bar_true / I_bar_true)
    }






  })

}

