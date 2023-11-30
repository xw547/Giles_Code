library("forestError")
# load data
data(airquality)
# remove observations with missing predictor variable values
airquality <- airquality[complete.cases(airquality), ]
# get number of observations and the response column index
n <- nrow(airquality)
response.col <- 1

# split data into training and test sets
train.ind <- sample(1:n, n * 0.9, replace = FALSE)
Xtrain <- airquality[train.ind, -response.col]
Ytrain <- airquality[train.ind, response.col]
Xtest <- airquality[-train.ind, -response.col]
# fit random forest to the training data
rf <- randomForest::randomForest(Xtrain, Ytrain, nodesize = 5,
                                 ntree = 500, keep.inbag = TRUE)
# compute out-of-bag prediction errors and locate each
# training observation in the trees for which it is out
# of bag
train_nodes <- findOOBErrors(rf, Xtrain)
# estimate conditional mean squared prediction errors,
# biases, prediction intervals, and error distribution
# functions for the test observations. provide
# train_nodes to avoid recomputing that step.
output <- quantForestError(rf, Xtrain, Xtest,
                           train_nodes = train_nodes)

output <- quantForestError(rf, Xtrain, Xtest,
                           what = c("p.error", "q.error"))

