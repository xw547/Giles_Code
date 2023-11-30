library("permimp")
library("party")
library("randomForest")
library("MASS")

### Data Generation
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
n = 1000

X = mvrnorm(n, rep(0, 12), Sigma)
y = X %*% beta + rnorm(n, 0, 0.5)
data_curr = data.frame(y, X)
X.test = mvrnorm(round(n/3), rep(0,12), Sigma)
y.test = X.test %*% beta + rnorm(round(n/3), 0, 0.5)
test_data = data.frame(y.test, X.test)



### Conducting Inference 
rf_original <- randomForest(y~., data = data_curr, keep.forest = T, keep.inbag = T)
# rf_original <- ctree(y~., data = data_curr)

CI <- permimp(rf_original, conditional = T, do_check = F)
VI <- permimp(rf_original, conditional = F, do_check = F)

mean((predict(rf_original,  test_data[,]) - test_data[,1])^2)
### With this package, I believe that the speed can be much faster than I had
### expected.


# Now, I'll try to calculate the EIF of conditional variable importance side. 



