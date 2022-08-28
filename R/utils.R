# ks gets the kernel estimation return yy.test=E[yy|xx=xx.test]
ks <- function(xx, yy, xx.test){
  nobs <- nrow(as.matrix(xx))
  nvars <- ncol(as.matrix(xx))
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(as.matrix(xx))==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    weighted.mean(yy, weight) # need to be changed mean(yy * weight) / (N-n)
  }
  if (is.null(dim(xx.test)) || NROW(xx.test)==1) {
    yy.test <- wm(xx.test)
  } else {
    if (nvars==1){
      yy.test <- sapply(as.matrix(xx.test), function(t){
        wm(t)
      })
    } else {
      yy.test <- apply(as.matrix(xx.test),1,function(t){
        wm(t)
      })
    }
  }
  yy.test
}

ks.missing <- function(xx, yy, samp_size, xx.test){ # rr should be T/F values
  nobs <- nrow(as.matrix(xx))
  nvars <- ncol(as.matrix(xx))
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(as.matrix(xx))==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    sum(yy * weight) / samp_size

  }
  if (is.null(dim(xx.test)) || NROW(xx.test)==1) {
    yy.test <- wm(xx.test)
  } else {
    if (nvars==1){
      yy.test <- sapply(as.matrix(xx.test), function(t){
        wm(t)
      })
    } else {
      yy.test <- apply(as.matrix(xx.test),1,function(t){
        wm(t)
      })
    }
  }
  yy.test
}


## Propensity score
fun.ps.model <- function(X, Z, R, model.type){
  if(sum(R==0)> 0){
    if(model.type=='glm'){
      glm.data <- data.frame(R=R, X=X, Z=Z)
      model <- glm(R ~ .-1, data = glm.data, family = "binomial")
    }else if(model.type=='glmnet'){
      cv.lasso.model <- cv.glmnet(x = cbind(X, Z=Z), y = R, family = "binomial", intercept=F)
      # cv.lasso.model
      coef(cv.lasso.model, s="lambda.min")
      model <- glmnet(x = cbind(X, Z= Z), y = R, family = "binomial", intercept = F, lambda=cv.lasso.model$lambda.min)
    }

    ## Conditional expectation with kernel smoother regression
  }else { model <- NULL }
  return(model)
}
get.ps <- function(model=NULL, X, Z){
  if(!is.null(model)){
    if(sum(class(model)=="glm")>0){
      ps <- predict(model, newdata = data.frame(X=X, Z=Z), type='response')

    }else if(sum(class(model)=="glmnet")>0){
      ps <- predict(model, newx = cbind(X, Z), type = "response")
    }else ps <- NULL
  }
  else ps <- NULL
  return(ps)
}

## Get nuisance model

getNuisance <- function(data, outcome=c('R', 'y'), method=c('lm', 'glmnet', 'kernel', 'density.ratio', 'true'), dimension.reduction = c("None", 'true', 'Screening', "slicedIR"),
                        sampleSplitIndex = NULL, Formula = NULL, predictAll = FALSE, screening.method="SIRS"){
  # if outcome = R, provide propensity score; if outcome = y, provide the conditional expectation
  data$predictor <- cbind(data$X, Z=data$Z)
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]

  # if sampleSplitIndex == null, samples are split
  if(is.null(sampleSplitIndex) & method %in% c("kernel", "lm", "density.ratio")){
    sampleSplitIndex <- rep(T, size)
    sampleSplitIndex[sample.int(size,0.5*size)] <- F
  }else if(is.null(sampleSplitIndex) & method %in% c("glmnet", "true")){
    sampleSplitIndex <- rep(T, size)
  }

  fit <- NULL
  supp <- rep(TRUE, times = p)
  supp.beta <- NULL
  dataTrain <- list(predictor=data$predictor[sampleSplitIndex,], outcome=data[[outcome]][sampleSplitIndex]) # Use sampleSplitIndex == TRUE for training
  dataPredict <- data$predictor[!sampleSplitIndex,] # Use sampleSplitIndex == FALSE or use all for predict
  if (predictAll){
    dataPredict=data$predictor
  }
  prediction <- NULL

  ### dimension.reduction ####
  # Screening using only training
  if(dimension.reduction=="Screening"){
    ans <- screening(dataTrain$predictor, dataTrain$outcome, method = screening.method, family = 'binomial')
    if (0.05*size >= p){
      supp <- (ans <= 5)
    }else{
      supp <- (ans <= floor(p/2))
    }
  } else if(dimension.reduction == "SlicedIR"){
    supp.beta <- LassoSIR::LassoSIR(X=dataTrain$predictor, Y=factor(dataTrain$outcome), categorical = T, screening = TRUE, no.dim=1)
    supp <- abs(supp.beta$beta)>0
  } else if(dimension.reduction == "true"){
    if(outcome == "R") {
      supp <- c(rep(c(T,F), 4), rep(F,p-8-1))
    }else if(outcome == "y"){
      supp <- c(rep(T,4), rep(F, p-4-1))
    }
  }


  # Building prediction model
  Formula <- function(support){
    expr <- (outcome ~ predictor)
    if (sum(support)==0){
      expr <- (outcome ~ 1)
    }
    expr
  }
  fit <- NULL

  dataTrain <- NULL
  if(dimension.reduction == "SlicedIR"){
    dataTrain=list(predictor = data$predictor[sampleSplitIndex,] %*% supp.beta$beta , outcome=data[[outcome]][sampleSplitIndex])
  }else{
    dataTrain=list(predictor = data$predictor[sampleSplitIndex,supp], outcome=data[[outcome]][sampleSplitIndex])
  }

  dataPredict <- NULL
  if(dimension.reduction == "SlicedIR"){
    dataPredict<-data$predictor %*% supp.beta$beta
  }else{
    dataPredict<-data$predictor[,supp]
  }
  if(!predictAll){
    dataPredict <- dataPredict[!sampleSplitIndex,]
  }




  ### method using the training set ####
  # available methods: lm, glmnet, kernel, density.ratio, true
  # kernel
  if ((method == 'lm')){
    if (sum(supp) > 0){
      fit <- glm(Formula(supp), family=binomial, data = dataTrain)
      prediction <- predict(fit, newdata = list(predictor=dataPredict), type="response")
    }
  } else if (method == "glmnet"){
    if (sum(supp) > 0){
      fit <- glmnet::cv.glmnet(x = dataTrain$predictor, y = dataTrain$outcome, family='binomial')
      prediction <- predict(fit, dataPredict, s='lambda.min',type='response' )
    }
  } else if (method == 'kernel') {
    if (sum(supp) > 0){
      prediction <- ks(dataTrain$predictor, dataTrain$outcome, dataPredict)
      prediction <- (prediction > 0.9) * 0.9 + (prediction < 0.1) * 0.1 + (prediction < 0.9) * (prediction > 0.1) * prediction
    }
  } else if (method == "density.ratio"){
    if (sum(supp) >0 ){
      f1 <- dataTrain$predictor[dataTrain$outcome==1,]
      f0 <- dataTrain$predictor[dataTrain$outcome==0,]
      w <- densratio::densratio(f1, f0, method = "uLSIF")
      rho <- mean(dataTrain$outcome)
      w_pred <- w$compute_density_ratio(dataPredict)

      prediction <- w_pred*rho/(w_pred*rho +(1-rho))
    }
  }

  prediction
}

screening <- function(x, y, method='glmnet', family='Gaussian'){
  var <- apply(x, 2, sd)
  supp <- order(var, decreasing = TRUE)
  if (method=='glmnet'){
    fit <- glmnet::cv.glmnet(x, y, family = family)
    coef <- fit$glmnet.fit$beta[,fit$lambda==fit$lambda.min]
    supp <- (abs(coef)>0)
  } else {
    fit <- VariableScreening::screenIID(x, y, method=method)
    supp <- fit$rank
  }
  supp
}

truncated <- function(values, thres = 0.1){
  if(sum(values<thres) >0 ) values[values < thres,] <- thres
  # if(sum(values>0.9) > 0) values[values > 0.9,] <- 0.9
  return(values)
}


loss_1st_boot <- function(samples, nuisance, beta_est, y_type, intercept){

  if(intercept){
    X <- cbind(1,samples$X)
  }else{
    X <- samples$X
    beta_est <- beta_est[-1,]
  }
  R <- samples$R; y <- samples$y;
  missingness <- nuisance$missingness; impute <- nuisance$impute;

  # t(X) %*% ((impute  - y )*R/missingness) + t(X) %*% ( b_1st(X, beta_est, y_type)-impute)
  t(X[R==1,]) %*% matrix((impute[R==1]  - y[R==1] )/missingness[R==1]) + t(X) %*% ( b_1st(X, beta_est, y_type)-impute)

}

# b_1st <- function(X, beta_tilde){
#
# }
b_1st <- function(X, beta_tilde, y_type){
  if(y_type=="binary"){
    1/(1+exp(-X %*% beta_tilde))
  }else if (y_type == "continuous"){
    X %*% beta_tilde
  }
}

b_2nd <- function(X, beta_tilde, y_type){
  if(y_type == "binary"){
    exp(-X %*% beta_tilde)/(1+exp(-X %*% beta_tilde))^2
  }else if(y_type == "continuous"){
    rep(1, nrow(X))
  }
}


fun_I <- function(X, interest_var, beta_tilde, w_hat, intercept, y_type){
  if(intercept){
    X <- cbind(1, X); interest_var <- interest_var + 1
    mean(as.matrix(b_2nd(X, beta_tilde, y_type) * X[,interest_var]* (X[,interest_var] -X[,-interest_var] %*% w_hat )))
  }else{
    mean(b_2nd(X, beta_tilde[-1,], y_type) * X[,interest_var]* (X[,interest_var] -X[,-interest_var] %*% w_hat ))
  }


}

make_sub <- function(samples, sampleSplitIndex, sample.div = T){
  samples_sub <- list()
  samples_sub$X <- samples$X[sampleSplitIndex == sample.div,] ;
  if(NCOL(samples$Z)>1){
    samples_sub$Z <- samples$Z[sampleSplitIndex == sample.div,];
  }else{
    samples_sub$Z <- samples$Z[sampleSplitIndex == sample.div];
  }

  samples_sub$R <- samples$R[sampleSplitIndex == sample.div];
  samples_sub$y <- samples$y[sampleSplitIndex == sample.div]
  samples_sub
}

