# R Package: missInfer

# Efficient surrogate-assisted estimation and inference in the presence of complex missing mechanisms 

Version: 1.0.0

Author: Jaeyoung Park <jaeyoungp@uchicago.edu>

Maintainer: Jaeyoung Park <jaeyoungp@uchicago.edu>

Description: The package implements a surrogate-assisted estimation in the presence of complex missing mechanisms. The package provides not only the propensity method but also our proposed method.

Imports: 
          glmnet, 
          mave
          
# Example
```
### Get details using the help function
?missInfer

### Generate data
set.seed(1823)
n <- 200; p <- 8; interest_lst <- 1:p; y_type <- "binary"; alpha <- -1;
thres_w <- 0.01; missing.rate <- 0.5
samples <- sample.generation(p = p, n= n, y_type = y_type, alpha = alpha, missing.rate =missing.rate)
X <- samples$X; Z <- samples$Z; y <- samples$y # R <- samples$R;

### Obtain debiased estimators
missInfer(X = X, Z = Z, y = y, interest_var = interest_lst, intercept = intercept, thres_w =  thres_w, pi.method = "kernel", dim.reduction = "onlyQ")
```

# Reference
Jaeyoung Park, Muxuan Liang, Yingqi Zhao, Xiang Zhong (2022). Efficient surrogate-assisted inference for patient-reported outcome with complex missing mechanism.
