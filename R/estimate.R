

fun_estimate <- function(X, Z=NULL, y, interest_var = 1:NCOL(X), intercept = T, thres_w = 0.01, variance = F){
  R <- (!is.na(y))*1
  p <- NCOL(X)

  conditioning <- ifelse(!is.null(Z), "ZX", "X")

  h_type <- ifelse(length(unique(y)== 2), "binary", "continuous")

  ZX <- cbind(Z, X)
  # mod_missingness <- glmnet::cv.glmnet(x = X, y = R, intercept = intercept, family="binomial")
  # supp.pi.glmnet <- glmnet::coef.glmnet(mod_missingness, s= 'lambda.min')
  #
  # dim.reduction.pi.X <- MAVE::mave(R ~ ., data = data.frame(R=R, X), method ="KSIR") #KSIR
  # cv.mave.pi.X <- MAVE::mave.dim(dim.reduction.pi.X)
  # supp.pi.X <- as.matrix(dim.reduction.pi.X$dir[[cv.mave.pi.X$dim.min]]) # True dimension

  dim.reduction.Q <- MAVE::mave(y ~ . , data = data.frame(ZX[R==1,], y = y[R==1]), method = "KSIR") #"KSIR"
  cv.mave <- MAVE::mave.dim(dim.reduction.Q)
  supp.Q <- as.matrix(dim.reduction.Q$dir[[cv.mave$dim.min]])

  samples <- list(n = NROW(X), p = p, X = as.matrix(X),
                     Y1 = Z,
                     Y2 = NA,  R = R, fun.h=y)

  imputed.Q <- predict_imputeQ(training = samples, test = samples,  #training = samples_tr
                                 conditioning = conditioning, supp.Q = supp.Q.X)

  beta_initial <- get_initial(samples = samples, impute = imputed.Q, intercept = intercept)

  # results_beta_x_sep_boot <- fun_estimate(samples = samples_tr, interest_lst = interest_lst,
  #                                         intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "X", pi.method = 'kernel', dim.reduction = "separate",
  #                                         supp.Q = supp.Q.X, supp.pi = supp.pi.X, beta_initial = beta_initial_X)
  # results_beta_x_sep_boot <- unlist(results_beta_x_sep_boot)
  # names(results_beta_x_sep_boot) <- paste0("x_sep_",names(results_beta_x_sep_boot), "_", rep(unique(c((!intercept)*1, interest_lst)), each=2))
  #
  # results_glmnet_boot <- fun_estimate(samples = samples_tr,interest_lst = interest_lst, #interest_lst,
  #                                     intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "X", pi.method = "glmnet", dim.reduction = "onlyQ",
  #                                     supp.Q = supp.Q.X, supp.pi = supp.pi.glmnet, beta_initial = beta_initial_X)
  # results_glmnet_boot <- unlist(results_glmnet_boot)
  # names(results_glmnet_boot) <- paste0("glmnet_",names(results_glmnet_boot), "_", rep(unique(c((!intercept)*1, interest_lst)), each = 2))

  beta_debiased <- fun_debias(samples = samples, interest_lst = interest_var,
                              intercept = intercept, thres_w =  thres_w, variance = F, conditioning = conditioning, pi.method = 'kernel', dim.reduction = "onlyQ",
                              supp.Q = supp.Q, supp.pi = supp.Q, beta_initial = beta_initial) #
  beta_debiased <- unlist(beta_debiased)
  names(beta_debiased) <- paste0("x_",names(beta_debiased), "_", rep(unique(c((!intercept)*1, interest_lst)), each=2))

  # results_beta_xy_boot <- fun_estimate(samples = samples_tr,  interest_lst = interest_lst,
  #                                      intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "Y1X", pi.method = 'kernel', dim.reduction = 'onlyQ',
  #                                      supp.Q = supp.Q.Y1X, supp.pi = supp.Q.Y1X, beta_initial = beta_initial_Y1X) #
  # results_beta_xy_boot <- unlist(results_beta_xy_boot)
  # names(results_beta_xy_boot) <- paste0("xy_",names(results_beta_xy_boot), "_", rep(unique(c((!intercept)*1, interest_lst)), each=2))

  if(intercept){
    c(beta_initial=as.numeric(beta_initial),
      beta_debiased[paste0("x_estimated_nuisance_",0:p)],
      # dim.pi = NCOL(supp.pi), dim.pi.glmnet = sum(supp.pi.glmnet!=0),
      dim.red = NCOL(supp.Q))
  }else{
    c(beta_initial=as.numeric(beta_initial[-1,]),
      beta_debiased[paste0("x_estimated_nuisance_",1:p)],
      # dim.pi = NCOL(supp.pi.X), dim.pi.glmnet = sum(supp.pi.glmnet!=0),
      dim.red.X = NCOL(supp.Q))
  }

}




set.seed(index)
idx_tr <- sample(1:nrow(demo_PRO9_5), size = round(0.5*nrow(demo_PRO9_5)))
idx_te <- setdiff(1:nrow(demo_PRO9_5), idx_tr)

tr <- demo_PRO9_5[idx_tr,]
te <- demo_PRO9_5[idx_te,]



# X_tr <- as.matrix(tr[, -which(colnames(tr) %in% c("MCID", "LOS.cat1", "LOS.cat2"))])
# X_te <- as.matrix(te[, -which(colnames(te) %in% c("MCID", "LOS.cat1", "LOS.cat2"))])
X_tr <- as.matrix(tr[, -which(colnames(tr) %in% c("MCID", "MCID.sur", "LOS.cat1", "LOS.cat2", "combined_sur", "No.cancel.noshow.post.cat1", "No.follow.up.cat1", "readmission1", "PRE_SURG_PROMIS_GLOBAL_HEALTH_T_SCORE_PHYSICAL",
                                                  'survey1.1', "survey1.2", "survey2.1", "survey2.2", "survey3.1", "survey3.2"))])
X_te <- as.matrix(te[, -which(colnames(te) %in% c("MCID", "MCID.sur","LOS.cat1", "LOS.cat2", "combined_sur", "No.cancel.noshow.post.cat1", "No.follow.up.cat1", "readmission1", "PRE_SURG_PROMIS_GLOBAL_HEALTH_T_SCORE_PHYSICAL",
                                                  'survey1.1', "survey1.2", "survey2.1", "survey2.2", "survey3.1", "survey3.2"))])
# X_tr <- as.matrix(tr[, -which(colnames(tr) %in% c("MCID", "MCID.sur", "No.cancel.noshow.post.cat", "No.cancel.noshow.post"))])
# X_te <- as.matrix(te[, -which(colnames(te) %in% c("MCID", "MCID.sur","No.cancel.noshow.post.cat", "No.cancel.noshow.post"))])
# X_tr <- as.matrix(tr[, -which(colnames(tr) %in% c("MCID", "MCID.sur", "No.follow.up"))])
# X_te <- as.matrix(te[, -which(colnames(te) %in% c("MCID", "MCID.sur","No.follow.up"))])
# X_tr <- as.matrix(tr[, -which(colnames(tr) %in% c("MCID", "No.cancel.noshow.post.cat"))])
# X_te <- as.matrix(te[, -which(colnames(te) %in% c("MCID", "No.cancel.noshow.post.cat"))])

R_tr <- (!is.na(tr$MCID))*1; R_te <- (!is.na(te$MCID))*1

p <- NCOL(X_tr)
intercept <- TRUE; thres_w <- 0.2; interest_lst <- 1:p; # intercept <- FALSE
h_type <- "indicator"; # h_type <- "difference"

mod_missingness <- cv.glmnet(x = X_tr, y = R_tr, intercept = intercept, family="binomial")
supp.pi.glmnet <- coef(mod_missingness, s= 'lambda.min')

dim.reduction.pi.X <- mave(R ~ ., data = data.frame(R=R_tr, X_tr), method ="KSIR") #KSIR
cv.mave.pi.X <- mave.dim(dim.reduction.pi.X)
supp.pi.X <- as.matrix(dim.reduction.pi.X$dir[[cv.mave.pi.X$dim.min]]) # True dimension

# dim.reduction.Q.X <- mave(MCID ~ . - LOS.cat2, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - LOS.cat1 - LOS.cat2, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - LOS.cat2 - No.cancel.noshow.post.cat1 - No.follow.up.cat1 , data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - PRE_SURG_PROMIS_GLOBAL_HEALTH_T_SCORE_PHYSICAL , data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - LOS.cat2 - No.cancel.noshow.post.cat1 - No.follow.up.cat1 - readmission1, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
dim.reduction.Q.X <- mave(MCID ~ . - combined_sur, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - survey1.1 - survey1.2 - survey2.1 - survey2.2 - survey3.1 - survey3.2, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - survey1.1 - survey1.2, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID.sur ~ . - MCID - LOS.cat1 - LOS.cat2, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID.sur ~ . - MCID - No.cancel.noshow.post.cat, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - No.cancel.noshow.post.cat, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - No.cancel.noshow.post, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.X <- mave(MCID ~ . - No.follow.up, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
cv.mave.X <- mave.dim(dim.reduction.Q.X)
supp.Q.X <- as.matrix(dim.reduction.Q.X$dir[[cv.mave.X$dim.min]])
# supp.Q.X <- as.matrix(dim.reduction.Q.X$dir[[10]])

dim.reduction.Q.Y1X <- mave(MCID ~ ., data = tr[R_tr==1,], method = "KSIR") #"KSIR"
# dim.reduction.Q.Y1X <- mave(MCID.sur ~ . - MCID, data = tr[R_tr==1,], method = "KSIR") #"KSIR"
cv.mave.Y1X <- mave.dim(dim.reduction.Q.Y1X)
supp.Q.Y1X <- as.matrix(dim.reduction.Q.Y1X$dir[[cv.mave.Y1X$dim.min]])
# supp.Q.Y1X <- as.matrix(dim.reduction.Q.Y1X$dir[[10]])




### performance ###
# setwd("/Users/jaeyoungpark/OneDrive - University of Florida/TJR readmissions/RcppEigen/")
source("R/Estimate_real_world.R")
source("R/utils.R")
# library(mgcv); library(glmnet); library(parallel); library(MAVE)


samples_tr <- list(n = nrow(X_tr), p = ncol(X_tr), X = as.matrix(X_tr),
                   # Y1 = tr[, c("No.follow.up")],
                   # Y1 = tr[, c("No.cancel.noshow.post")],
                   # Y1 = tr[, c("No.cancel.noshow.post.cat")],
                   # Y1 = tr[, c("LOS.cat1", "LOS.cat2")],
                   # Y1 = tr[, c("LOS.cat2")],
                   # Y1 = tr[, c("LOS.cat2", "No.follow.up.cat1", "No.cancel.noshow.post.cat1", "readmission1")],
                   # Y1 = tr[, c("LOS.cat2", "No.follow.up.cat1", "No.cancel.noshow.post.cat1")],
                   # Y1 = tr[, c("PRE_SURG_PROMIS_GLOBAL_HEALTH_T_SCORE_PHYSICAL")],
                   # Y1 = tr[, c('survey1.1', "survey1.2", "survey2.1", "survey2.2", "survey3.1", "survey3.2")],
                   # Y1 = tr[, c('survey1.1', "survey1.2")],
                   Y1 = tr[, c("combined_sur")],
                   Y2 = NA,  R = R_tr, fun.h=tr$MCID)


imputed.Q.X <- predict_imputeQ(training = samples_tr, test = samples_tr,  #training = samples_tr
                               conditioning = "X", supp.Q = supp.Q.X)

imputed.Q.Y1X <- predict_imputeQ(training = samples_tr, test = samples_tr,#training = samples_tr
                                 conditioning = "Y1X", supp.Q = supp.Q.Y1X)

beta_initial_X <- get_initial(samples = samples_tr, impute = imputed.Q.X, intercept = intercept)
beta_initial_Y1X <- get_initial(samples = samples_tr, impute = imputed.Q.Y1X, intercept = intercept)



results_beta_x_sep_boot <- fun_estimate(samples = samples_tr, interest_lst = interest_lst,
                                        intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "X", pi.method = 'kernel', dim.reduction = "separate",
                                        supp.Q = supp.Q.X, supp.pi = supp.pi.X, beta_initial = beta_initial_X)
results_beta_x_sep_boot <- unlist(results_beta_x_sep_boot)
names(results_beta_x_sep_boot) <- paste0("x_sep_",names(results_beta_x_sep_boot), "_", rep(unique(c((!intercept)*1, interest_lst)), each=2))

results_glmnet_boot <- fun_estimate(samples = samples_tr,interest_lst = interest_lst, #interest_lst,
                                    intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "X", pi.method = "glmnet", dim.reduction = "onlyQ",
                                    supp.Q = supp.Q.X, supp.pi = supp.pi.glmnet, beta_initial = beta_initial_X)
results_glmnet_boot <- unlist(results_glmnet_boot)
names(results_glmnet_boot) <- paste0("glmnet_",names(results_glmnet_boot), "_", rep(unique(c((!intercept)*1, interest_lst)), each = 2))


# source("R/Estimate_112221_revised.R")
# set.seed(1)
results_beta_x_boot <- fun_estimate(samples = samples_tr, interest_lst = interest_lst,
                                    intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "X", pi.method = 'kernel', dim.reduction = "onlyQ",
                                    supp.Q = supp.Q.X, supp.pi = supp.Q.X, beta_initial = beta_initial_X) #
results_beta_x_boot <- unlist(results_beta_x_boot)
names(results_beta_x_boot) <- paste0("x_",names(results_beta_x_boot), "_", rep(unique(c((!intercept)*1, interest_lst)), each=2))

results_beta_xy_boot <- fun_estimate(samples = samples_tr,  interest_lst = interest_lst,
                                     intercept = intercept, thres_w =  thres_w, variance = F, conditioning = "Y1X", pi.method = 'kernel', dim.reduction = 'onlyQ',
                                     supp.Q = supp.Q.Y1X, supp.pi = supp.Q.Y1X, beta_initial = beta_initial_Y1X) #
results_beta_xy_boot <- unlist(results_beta_xy_boot)
names(results_beta_xy_boot) <- paste0("xy_",names(results_beta_xy_boot), "_", rep(unique(c((!intercept)*1, interest_lst)), each=2))

no_debiased <- cv.glmnet(x = X_tr[R_tr ==1,], y = tr$MCID[R_tr ==1], family='binomial', intercept = intercept)


cbind(results_beta_x_sep_boot[paste0("x_sep_estimated_nuisance_",unique(c((!intercept)*1, 1:p)))],
      results_glmnet_boot[paste0("glmnet_estimated_nuisance_",unique(c((!intercept)*1, 1:p)))],
      results_beta_x_boot[paste0("x_estimated_nuisance_",unique(c((!intercept)*1, 1:p)))],
      results_beta_xy_boot[paste0("xy_estimated_nuisance_",unique(c((!intercept)*1, 1:p)))])
# coef(no_debiased, s = "lambda.min")[ifelse(intercept, 1, 2):(p+1)])


samples_te <- list(n = nrow(X_te), p = ncol(X_te), X = as.matrix(X_te),
                   # Y1 = te[, c("No.follow.up")],
                   # Y1 = te[, c("No.cancel.noshow.post")],
                   # Y1 = te[, c("No.cancel.noshow.post.cat")],
                   # Y1 = te[, c("LOS.cat1", "LOS.cat2")],
                   # Y1 = te[, c("LOS.cat2")],
                   Y1 = te[, c("combined_sur")],
                   # Y1 = te[, c('survey1.1', "survey1.2", "survey2.1", "survey2.2", "survey3.1", "survey3.2")],
                   # Y1 = te[, c("LOS.cat2", "No.follow.up.cat1", "No.cancel.noshow.post.cat1", "readmission1")],
                   # Y1 = te[, c("LOS.cat2", "No.follow.up.cat1", "No.cancel.noshow.post.cat1")],
                   # Y1 = te[, c("PRE_SURG_PROMIS_GLOBAL_HEALTH_T_SCORE_PHYSICAL")],
                   # Y1 = te[, c('survey1.1', "survey1.2")],
                   Y2 = NA,  R = R_te, fun.h=te$MCID)



### Need to change??

pred.Q.X <- predict_imputeQ(training = samples_te, test = samples_te,  #training = samples_tr
                            conditioning = "X",
                            supp.Q = supp.Q.X)

pred.Q.Y1X <- predict_imputeQ(training = samples_te, test = samples_te,#training = samples_tr
                              conditioning = "Y1X",
                              supp.Q = supp.Q.Y1X)
summary(cbind(pred.Q.X, pred.Q.Y1X))

# pred.funh1 <- X_te %*% results_beta_x_sep_boot[paste0("x_sep_estimated_nuisance_",1:p)]
# pred.funh2 <- X_te %*% results_glmnet_boot[paste0("glmnet_estimated_nuisance_",1:p)]
# pred.funh3 <- X_te %*% results_beta_x_boot[paste0("x_estimated_nuisance_",1:p)]
# pred.funh4 <- X_te %*% results_beta_xy_boot[paste0("xy_estimated_nuisance_",1:p)]
# summary(cbind(pred.funh1, pred.funh2, pred.funh3, pred.funh4))



if(intercept){
  # pred.funh0 <- cbind(1,X_te) %*% coef(no_debiased, s = "lambda.min")[1:(p+1)]
  pred.funh1 <- cbind(1,X_te) %*% results_beta_x_sep_boot[paste0("x_sep_estimated_nuisance_",0:p)]
  pred.funh2 <- cbind(1,X_te) %*% results_glmnet_boot[paste0("glmnet_estimated_nuisance_",0:p)]
  pred.funh3 <- cbind(1,X_te) %*% results_beta_x_boot[paste0("x_estimated_nuisance_",0:p)]
  pred.funh4 <- cbind(1,X_te) %*% results_beta_xy_boot[paste0("xy_estimated_nuisance_",0:p)]
}else{
  # pred.funh0 <- X_te %*% coef(no_debiased, s = "lambda.min")[2:(p+1)]
  pred.funh1 <- X_te %*% results_beta_x_sep_boot[paste0("x_sep_estimated_nuisance_",1:p)]
  pred.funh2 <- X_te %*% results_glmnet_boot[paste0("glmnet_estimated_nuisance_",1:p)]
  pred.funh3 <- X_te %*% results_beta_x_boot[paste0("x_estimated_nuisance_",1:p)]
  pred.funh4 <- X_te %*% results_beta_xy_boot[paste0("xy_estimated_nuisance_",1:p)]
}

summary(cbind( pred.funh1, pred.funh2, pred.funh3, pred.funh4))

if(h_type == "difference"){
  res.dev <- c(
    # b0 = 1/test.n *(0.5 * t(pred.funh0 ) %*% pred.funh0 - t(test.set$fun.h) %*% pred.funh0),
    b1 = mean(0.5 * pred.funh1  * pred.funh1 - pred.Q.X * pred.funh1),
    b2 = mean(0.5 * pred.funh2  * pred.funh2 - pred.Q.X * pred.funh2),
    p1 = mean(0.5 * pred.funh3  * pred.funh3 - pred.Q.X * pred.funh3),
    p2 = mean(0.5 * pred.funh4  * pred.funh4 - pred.Q.Y1X * pred.funh4))

  # res.dev <- c(
  #   # b0 = 1/test.n *(0.5 * t(pred.funh0 ) %*% pred.funh0 - t(test.set$fun.h) %*% pred.funh0),
  #   b1 = mean(0.5 * pred.funh1  * pred.funh1 - te$MCID[R_te==1] * pred.funh1),
  #   b2 = mean(0.5 * pred.funh2  * pred.funh2 - te$MCID[R_te==1] * pred.funh2),
  #   p1 = mean(0.5 * pred.funh3  * pred.funh3 - te$MCID[R_te==1] * pred.funh3),
  #   p2 = mean(0.5 * pred.funh4  * pred.funh4 - te$MCID[R_te==1] * pred.funh4))

}else if(h_type == "indicator"){
  res.dev <- c(
    # b0 = mean(log(1+exp(pred.funh0)) - test.set$fun.h * pred.funh0),
    b1 = mean(log(1+exp(pred.funh1)) - pred.Q.X * pred.funh1),
    b2 = mean(log(1+exp(pred.funh2)) - pred.Q.X * pred.funh2),
    p1 = mean(log(1+exp(pred.funh3)) - pred.Q.X * pred.funh3),
    p2 = mean(log(1+exp(pred.funh4)) - pred.Q.Y1X * pred.funh4)
  )

  # res.dev <- c(
  #   # b0 = mean(log(1+exp(pred.funh0)) - test.set$fun.h * pred.funh0),
  #   b1 = mean(log(1+exp(pred.funh1)) - te$MCID[R_te==1] * pred.funh1),
  #   b2 = mean(log(1+exp(pred.funh2)) - te$MCID[R_te==1] * pred.funh2),
  #   p1 = mean(log(1+exp(pred.funh3)) - te$MCID[R_te==1] * pred.funh3),
  #   p2 = mean(log(1+exp(pred.funh4)) - te$MCID[R_te==1] * pred.funh4)
  # )
}


if(intercept){
  c(beta_initial=as.numeric(beta_initial_X),
    results_beta_x_sep_boot[paste0("x_sep_estimated_nuisance_",0:p)], results_glmnet_boot[paste0("glmnet_estimated_nuisance_",0:p)],
    results_beta_x_boot[paste0("x_estimated_nuisance_",0:p)], results_beta_xy_boot[paste0("xy_estimated_nuisance_",0:p)],
    dim.pi = NCOL(supp.pi.X), dim.pi.glmnet = sum(supp.pi.glmnet!=0), dim.red.X = NCOL(supp.Q.X), dim.red.Y1X = NCOL(supp.Q.Y1X),
    res.dev)
}else{
  c(beta_initial=as.numeric(beta_initial_X[-1,]),
    results_beta_x_sep_boot[paste0("x_sep_estimated_nuisance_",1:p)], results_glmnet_boot[paste0("glmnet_estimated_nuisance_",1:p)],
    results_beta_x_boot[paste0("x_estimated_nuisance_",1:p)], results_beta_xy_boot[paste0("xy_estimated_nuisance_",1:p)],
    dim.pi = NCOL(supp.pi.X), dim.pi.glmnet = sum(supp.pi.glmnet!=0), dim.red.X = NCOL(supp.Q.X), dim.red.Y1X = NCOL(supp.Q.Y1X),
    res.dev)
}

# })
}))
