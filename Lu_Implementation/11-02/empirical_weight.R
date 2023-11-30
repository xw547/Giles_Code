### This function returns a matrix of index and their corresponding 
### empirical probability.

### Notice that this method is designed to estimate the corresponding 
### conditional density of the samples IN the training set, we do need to
### modify the function, even the methodology when we are considering results
### outside the training set.
library(data.table)

empirical_rf_pdf <- function(forest, X.train, X.test){
  
  # get terminal nodes of training observations
  train.terminal.nodes <- attr(predict(forest, X.train, nodes = TRUE), "nodes")
  
  # get number of times each training observation appears in each tree
  bag.count <- forest$inbag
  
  # get the terminal nodes of the training observations in the trees in which
  # they are OOB (for all other trees, set the terminal node to be NA)
  train.terminal.nodes[bag.count != 0] <- NA
  
  # reshape train.terminal.nodes to be a long data.table and include OOB
  # prediction errors as a column
  train_nodes_emp <- data.table::melt(
    data.table::as.data.table(train.terminal.nodes)[, `:=`(obsid_train = .I)],
    
    measure.vars = 1:ncol(train.terminal.nodes),
    variable.name = "tree",
    value.name = c("terminal_node"),
    variable.factor = FALSE,
    na.rm = TRUE)
  
  # collapse the long data.table by unique tree/node
  train_nodes_emp <- train_nodes_emp[,
                             .(obsid_train = list(obsid_train)),
                             keyby = c("tree", "terminal_node")]
  
  # get test predictions
  test.preds <- predict(forest, X.test, nodes = TRUE)
  
  # get terminal nodes of test observations
  test.terminal.nodes <- attr(test.preds, "nodes")
  
  # format test observation predictions
  attr(test.preds, "nodes") <- NULL
  
  # reshape test.terminal.nodes to be a long data.table and
  # add unique IDs and predicted values
  test_nodes_emp <- data.table::melt(
    data.table::as.data.table(test.terminal.nodes)[, `:=`(testid_test = .I, pred = test.preds)],
    id.vars = c("testid_test", "pred"),
    measure.vars = 1:ncol(test.terminal.nodes),
    variable.name = "tree",
    variable.factor = FALSE,
    value.name = "terminal_node")
  
  # set key columns for faster indexing
  data.table::setkey(test_nodes_emp, tree, terminal_node)
  data.table::setkey(train_nodes_emp, tree, terminal_node)
  list_combined <- train_nodes_emp[test_nodes_emp,
              .(tree, terminal_node, testid_test, pred, obsid_train)]
  result_matrix = matrix(0, nrow = dim(X.test)[1], ncol = dim(X.train)[1])
  for (i in 1:dim(X.test)[1]){
    #print(i)
    entries = table_to_density(list_combined, i)
    #print(head(entries))
    result_matrix[i,as.numeric(names(entries))] = entries
  }
  return(result_matrix)
}

table_to_density <- function(list_combined, i){
  temp_list   = list_combined[testid_test == i]
  temp_result = unlist(temp_list$obsid_train)
  return(table(temp_result)/length(temp_result))
}


