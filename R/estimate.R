

#' Estimate coefficients under complex missing mechanisms
#'
#' @param X Input data matrix
#' @param Z Surrogate outcome(s) Z can be a vector or a matrix (default: NULL)
#' @param y Response variable
#' @param interest_var Variables to be debiased (default: the first four variables)
#' @param intercept Intercept (default: TRUE)
#' @param thres_w Threshold to avoid too small inverse pi value (default: 0.01)
#' @param pi.method Method for estimating the inverse pi function. The kernel method is default but the glmnet method is also available. Should be either "kernel" or "glmnet".
#' @param dim.reduction Models for estimating a reduced dimension. The imputation model is default but both nuisance models are used for the dimension reduction. Should be either "separate" or "onlyQ".
#'
#' @return a list that contains the initial coefficients, debaised coefficients, dimensionality of reduced dimension
#' @export
#'
#' @examples
#' # Generate data
#' set.seed(1823)
#' n <- 200; p <- 8; interest_lst <- 1:8; y_type <- "binary"; alpha <- -1;
#' thres_w <- 0.01; missing.rate <- 0.5
#' samples <- sample.generation(p = p, n= n, y_type = y_type, alpha = alpha, missing.rate =missing.rate)
#' X <- samples$X; Z <- samples$Z; y <- samples$y # R <- samples$R;
#'
#' # Obtain debiased estimators
#' missInfer(X = X, Z = Z, y = y, interest_var = interested_lst, intercept = intercept, thres_w =  thres_w, pi.method = "kernel", dim.reduction = "onlyQ")

missInfer <- function(X, Z=NULL, y, interest_var = 1:4, intercept = T, thres_w = 0.01,
                         pi.method = c("kernel", "glmnet"), dim.reduction = c("separate", "onlyQ")
                         ){
  if(sum(class(X) == "matrix") == 0 ) stop("X should be a matrix.")
  if(!(pi.method %in% c("kernel", "glmnet")) ) stop("pi.method should be either kernel or glmnet.")
  if(!(dim.reduction %in% c("separate", "onlyQ"))) stop("dim.reduction should be either separate or onlyQ.")
  if(pi.method == "glmnet") {
    dim.reduction <- "onlyQ"
    message("When pi.method is glmnet, the dimension reduction is conducted separately for pi and Q")
  }


  R <- (!is.na(y))*1
  p <- NCOL(X); n <- NROW(X)

  conditioning <- ifelse(!is.null(Z), "ZX", "X")

  y_type <- ifelse(length(unique(y)== 2), "binary", "continuous")

  ZX <- cbind(Z, X)


  dim.reduction.Q <- MAVE::mave(y ~ . , data = data.frame(ZX[R==1,], y = y[R==1]), method = "KSIR") #"KSIR"
  cv.mave <- MAVE::mave.dim(dim.reduction.Q)
  supp.Q <- as.matrix(dim.reduction.Q$dir[[cv.mave$dim.min]])

  # Specify supp.pi by the pi.method and the dim.reduction.
  supp.pi <- supp.Q #
  if(pi.method == "glmnet"){
    mod_missingness <- glmnet::cv.glmnet(x = ZX, y = R, intercept = intercept, family="binomial")
    supp.pi <- glmnet::coef.glmnet(mod_missingness, s= 'lambda.min')
  }
  if((pi.method == "kernel") & (dim.reduction == "separate")){
    dim.reduction.pi <- MAVE::mave(R ~ ., data = data.frame(R=R, ZX), method ="KSIR") #KSIR
    cv.mave.pi <- MAVE::mave.dim(dim.reduction.pi)
    supp.pi <- as.matrix(dim.reduction.pi$dir[[cv.mave.pi$dim.min]])
  }

  samples <- list(n = n, p = p, X = as.matrix(X),
                     Z = Z, R = R, y=y)

  imputed.Q <- predict_imputeQ(training = samples, test = samples,  #training = samples_tr
                                 conditioning = conditioning, supp.Q = supp.Q)

  beta_initial <- get_initial(samples = samples, impute = imputed.Q, intercept = intercept)

  beta_debiased <- fun_debias(samples = samples, interest_var = interest_var,
                              intercept = intercept, thres_w =  thres_w, variance = F, conditioning = conditioning, pi.method = pi.method, dim.reduction = dim.reduction,
                              supp.Q = supp.Q, supp.pi = supp.pi, beta_initial = beta_initial) #
  beta_debiased <- unlist(beta_debiased)
  names(beta_debiased) <- paste0(names(beta_debiased), "_", rep(unique(c((!intercept)*1, interest_var)), each=1))

  if(intercept){
    c(beta_initial=as.numeric(beta_initial),
      beta_debiased[paste0("debiased_",0:p)],
      dim.red = NCOL(supp.Q))
  }else{
    c(beta_initial=as.numeric(beta_initial[-1,]),
      beta_debiased[paste0("debiased_",1:p)],
      dim.red = NCOL(supp.Q))
  }

}

### get initial betas ####

get_initial <- function(samples, impute, intercept = F){

  X <- samples$X; Z <- samples$Z; R <- samples$R; y <- samples$y
  n <- nrow(X); p <- ncol(X);
  if(length(unique(y)[!is.na(unique(y))])==2 ){
    y_type <- "binary"
  }else { y_type <- "continuous"}


  count <- 1
  if(y_type == "continuous"){
    # mod_impute <- glmnet(x = X, y = impute_entire, intercept=intercept, lambda = 0)
    mod_impute <- glmnet::cv.glmnet(x = X, y = impute, intercept=intercept)
  }else if(y_type == "binary"){
    err_val <- 0; count <- 0
    while(err_val == 0){
      tryCatch({
        mod_impute <- glmnet::cv.glmnet(x = rbind(X, X), y = factor(c(rep(0,n), rep(1, n))), weights = c(1-impute, impute), intercept=intercept, family="binomial")
        err_val <- 1
      }, error = function(err){
        err_val <- 0;
        if(count >= 50){
          message("Resampling is needed. Here is the original error message:")
          message(err)
          return(NA)
        }
        count <- count + 1
      })
    }
  }
  beta_initial <- coef(mod_impute, s = "lambda.min")
  return(beta_initial)
}


#### Estimate nuisance functions #####

predict_imputeQ <- function(training, test, nuisance_true = NULL, conditioning = c("X", "ZX"),
                            supp.Q = NULL){

  #### variable assignments #####
  X <- training$X; Z <- training$Z; R <- training$R; y <- training$y
  n <- nrow(X); p <- ncol(X);

  X_te <- test$X; Z_te <- test$Z; R_te <- test$R; y_te <- test$y
  n_te <- nrow(X_te);


  #### model for missingness ####

  if(conditioning == "X"){
    combined_dataset <- X
    combined_test <- X_te
  }else if(conditioning == "ZX"){
    combined_dataset <- cbind(Z= Z, X)
    combined_test <- cbind(Z= Z_te, X_te)
  }



  ### Q function ####
  if(!is.null(supp.Q)){
    impute_Q <- ks(xx=as.matrix(combined_dataset[R==1,]) %*% supp.Q, yy = y[R==1], xx.test = as.matrix(combined_test) %*% supp.Q)

  }else{
    impute_Q <- ks(xx=as.matrix(combined_dataset[R==1,]), yy = y[R==1], xx.test = as.matrix(combined_tes))
  }

  return(impute_Q)
}


