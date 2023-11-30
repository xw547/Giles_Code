library(MASS)
library(ranger)

n = 1000

beta = c(5, 5, 2, 0, -5, -5, -2, 0, rep(0, 4))
# Initialize an empty matrix
Sigma <- matrix(0, nrow = 12, ncol = 12)

# Fill in the values where the absolute difference in row and column indices is <= 1
for (i in 1:4) {
  for (j in 1:12) {
    if (abs(i - j) <= 4) {
      Sigma[i, j] <- 0.9
    }
  }
}


diag(Sigma) <- rep(1,12)


### Data generation
set.seed(2023)

X = mvrnorm(n, rep(0, 12), Sigma)
y = X %*% beta + rnorm(n, 0, 0.5)
data_curr = data.frame(y, X)
X.test = mvrnorm(round(n/3), rep(0,12), Sigma)
y.test = X.test %*% beta + rnorm(round(n/3), 0, 0.5)
test_data = data.frame(y.test, X.test)



rf_origin  <- ranger(y~., data = data_curr)


### In the following part, I'll try to 
### 1. Calculate the permutation importance for one of the features
### 2. Calculate the conditional permutation importance 

##### Variable importance for some of the features: illustrative examples
##### Notice that this isn't quite how it should be played, since we should 
##### compare the performance of each tree instead.

X_1pi = randomize_column(X, 1)
rf_1pi <- ranger(y~., data = data.frame(y, X_1pi))

VI_1 = mean((predict(rf_origin, data = test_data[,-1])$predictions - test_data[,1])^2 -
  (predict(rf_1pi, data = test_data[,-1])$predictions - test_data[,1])^2)

X_2pi = randomize_column(X, 2)
rf_2pi <- ranger(y~., data = data.frame(y, X_2pi))

VI_2 = mean((predict(rf_origin, data = test_data[,-1])$predictions - test_data[,1])^2 -
              (predict(rf_2pi, data = test_data[,-1])$predictions - test_data[,1])^2)

  
X_10pi = randomize_column(X, 10)
rf_10pi <- ranger(y~., data = data.frame(y, X_10pi))

VI_10 = mean((predict(rf_origin, data = test_data[,-1])$predictions - test_data[,1])^2 -
              (predict(rf_10pi, data = test_data[,-1])$predictions - test_data[,1])^2)

##### Conditional variable importance


### Auxiliary results.
randomize_column <- function(matrix, col_index) {
    matrix[, col_index] <- sample(matrix[, col_index])
    return(matrix)
}






