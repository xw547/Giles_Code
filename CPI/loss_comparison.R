### The purpose of this script is to test how large the difference can be 
### as derived in the cpi part. 

### Loading required packages
library("MASS")
library("VGAM")
library("permimp")
library("party")
library("randomForest")
library("MASS")

### Problem Set-up
### I'll largely follow the formulation of Hooker and Zhou 2021 SC 

set.seed(2021)
n    = 1000
rho  = 0.5
X_12 = rbinormcop(n, rho)
X_310 = matrix(runif(n*8, 0,1), ncol = 8)
X = cbind(X_12, X_310)
beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)

y = X%*%beta + rnorm(n, 0, 1)

### Prediction and corresponding loss
rf_original <- randomForest(y~., data = data.frame(y,X) , keep.forest = T, 
                            keep.inbag = T)

length_ana_result <- length_analysis(rf_original, conditional = T)
ci_rf_original <- permimp(rf_original)


test_tree = getTree(rf_original)
getOOB(rf_original,1)
















