## permimp for randomForest: permimp.randomForest.R
loss.randomForest <- function (object, nperm = 1, OOB = TRUE, scaled = FALSE,
                               conditional = FALSE, threshold = .95, whichxnames = NULL,   
                               thresholdDiagnostics = FALSE, progressBar = TRUE, 
                               do_check = TRUE, ...)
{
  # get randomForest call
  rfCall <- match.call(randomForest, object$call, expand.dots = TRUE)
  
  # check if formula was used
  if(inherits(object, "randomForest.formula")){
    
    # select input and responses (y) - using the randomForest.fromula-call
    dfCall <- rfCall
    dfCall[[1L]] <- quote(stats::model.frame)
    dfCall <- dfCall[c(1, match(c("formula", "data", "subset", "weights", "na.action", 
                                  "offset"), names(dfCall), 0L))]
    
    # make a model frame
    mf <- eval(dfCall, parent.frame())   # model frame  
    y <- model.response(mf)       # extract response
    mt <- attr(mf, "terms")              # model terms
    attr(mt, "intercept") <- 0           # intercept should not be included
    input <- model.frame(terms(reformulate(attributes(mt)$term.labels)), data.frame(mf))
    rm("mf", "mt")
    
  } 
  else {
    
    # select input and responses (y) - using the randomForest-call
    y <- eval(rfCall$y, parent.frame())       # extract response
    input <- eval(rfCall$x, parent.frame())   # extract input
    
  }
  # weights etc. are not included in the current computation
  if(object$type == "classification" & !is.null(rfCall$classwt)) 
    warning(sQuote("permimp"), " does not take ", sQuote("classwt"), " into account during permuation and prediction. \n", "All observations are automatically weighed equally.",
            immediate. = TRUE)
  
  if(object$type == "classification" & !is.null(rfCall$cutoff)) 
    warning(sQuote("permimp"), " does not take ", sQuote("cutoff"), " into account during prediction. \n", "The default cut-off is automatically applied.",
            immediate. = TRUE)
  
  if(!is.null(rfCall$sampsize) & length(eval(rfCall$sampsize, parent.frame())) > 1) 
    warning(sQuote("permimp"), " is based on the OOB values, using stratification to sample the IB values \n", "may have an undesired impact.",
            immediate. = TRUE)
  
  out <- dolength(object, input, inp = NULL, y, OOB, threshold, conditional, 
                   whichxnames, ntree = object$ntree, nperm, scaled,
                   progressBar, thresholdDiagnostics, 
                   w = NULL, AUC = FALSE, 
                   pre1.0_0 = TRUE, mincriterion = NULL, asParty = FALSE)
  
  return(out)
  
}

## doPermimp is the working horse of the permimp methods. 
## is called by all the permimp methods.

dolength <- function(object, input, inp, y, OOB, threshold, conditional, 
                      whichxnames, ntree, nperm, scaled,
                      progressBar, thresholdDiagnostics, 
                      w, AUC, pre1.0_0, mincriterion, asParty)
{
  # Check if conditional permutation importance is possible
  if (conditional) {
    if(!all(complete.cases(input)))
      stop("cannot compute variable importance measure with missing values")
    
    if (conditional && threshold == 1) {
      warning(sQuote("permimp"), 
              paste0(": Unable to permute conditionally. \n",  
                     "The chosen threshold is too high, no variables to condition on were selected. \n",
                     "Instead the unconditional permimp values are computed. "),
              call. = FALSE, immediate. = TRUE)
      doPermimpCall <- match.call()
      doPermimpCall$conditional <- FALSE
      return(eval(doPermimpCall))
    }
  }
  
  # select the predictors for which to compute the permutation importance
  xnames <- colnames(input)
  if(is.null(whichxnames)) {
    whichxnames <- xnames
    whichVarIDs <- seq_along(xnames)
  }
  else {
    whichVarIDs <- match(whichxnames, table = xnames)
    if(length(whichVarIDs) < 1){
      stop("Error: whichxnames is not a subset of the predictor variable names in the forest.")
    }
    whichVarIDs <- whichVarIDs[order(whichVarIDs)]
  } 
  
  # Check outcome and selet the relevant error- and pred-functions
  type <- getOutcomeType(object)
  error <- selectError(type, AUC)
  nullError <- selectNullError(type)
  pred <- selectPred(object, type, w, inp, y)
  
  # when asParty == TRUE, collect cond list first
  if (conditional && asParty) {
    cond_list <- create_cond_list(binnedVars = NULL, threshold,
                                  input, seq_along(xnames), asParty = TRUE)
  }
  
  ## list for several permutations
  ## this array is initialized with values 0 so that a tree that does not 
  ## contain the current variable adds importance 0 to its average importance
  perror <- array(0, dim = c(ntree, length(xnames), nperm), 
                  dimnames = list(NULL, xnames, NULL))
  
  ## this matrix will be used to give suggestions to de/increase the used threshold
  ## it is initialized with values NA.
  changeThres <- array(NA, dim = c(ntree, length(xnames), nperm), 
                       dimnames = list(NULL, xnames, NULL))
  
  # start progress bar
  if(progressBar) pBar <- txtProgressBar(min = 0, max = ntree, 
                                         style = 3, char = "|")
  length_mat = matrix(NA, ncol = length(xnames), nrow = ntree)
  partition_mat = c()
  # for all trees (treeNr) in the forest
  for (treeNr in seq_len(ntree)){
    tree <- getTree(object, treeNr)
    
    ## if OOB == TRUE use only oob observations, otherwise use all observations in learning sample
    if(OOB){oob <- getOOB(object, treeNr)} else {oob <- rep(TRUE, length(y))}
    
    ## prediction & error before permutation
    p <- pred(tree, inp, mincriterion, -1L, input)
    eoob <- error(p, oob, y)
    
    ## select variables that are used for splitting in the current tree
    varsInTree <- intersect(unique(varIDs(tree)), whichVarIDs)
    
    ## Only make the binned variables based on splitting points when conditional == TRUE
    if(conditional) {
      ## make list of variables, categorized/binned using the used splitting points
      binnedVars <- makeBinnedVars(varsInTree, tree, oob, input) 
      if (!asParty) {
        cond_list <- create_cond_list(binnedVars, threshold, input, varsInTree, asParty = FALSE)
        length_mat[treeNr, ] <- vector_lengths(cond_list)
      }
    }
    for(j in varsInTree){
      
      if (!conditional && !pre1.0_0){
        ## splitwise permutation only possible for RandomForest (party) object.
        p <- pred(tree, inp, mincriterion, as.integer(j))
      }
      else {
        if(conditional){
          changeThres[treeNr, j, per] <- 0  
          # if variable is in tree, NA -> 0
          ## only condition on variables that are in tree, 
          ## and that are associated with the current variable
          if(asParty) 
            varsToCondOn <- intersect(cond_list[[as.character(j)]], varsInTree)
          else varsToCondOn <- cond_list[[as.character(j)]]
          
          if(length(varsToCondOn) < 1){
            ## If there are no variables to condition on, conditionally permuting is impossible.
            ## -1 corresponds to a suggestion to decrease the used threshold
            changeThres[treeNr, j, per] <- -1  
            perm <- sample(which(oob))
          } else {
            samples <- samplel_perm(varID = j, varsToCondOn,
                                     binnedVars, oob, asParty)
            partition_mat = c(partition_mat, samples)
            
          }
        }
      }
      
      
    } ## end of for(j in varsInTree)

  } ## end of for (treeNr in 1:ntree)
  
  return(list(partition_mat, length_mat))
}

samplel_perm<- function(varID, varsToCondOn, binnedVars, oob, asParty) 
{
  ## same results as party
  if(asParty){
    CondPartitions <- interaction(binnedVars[as.character(varsToCondOn)], drop = TRUE, sep = "")
    levels(CondPartitions) <- 1:nlevels(CondPartitions)
    parts <- listParts(CondPartitions, oob)
  } 
  else {
    CondPartitions <- fastInteraction(binnedVars[as.character(varsToCondOn)])
    levels(CondPartitions) <- 1:nlevels(CondPartitions)
    ## not exactly the same results as party
    ## make partitions (including variable of interest)
    ## Not very sure why the all partitions should exist 
    AllPartitions <- fastInteraction(list(CondPartitions, binnedVars[names(binnedVars) == varID][[1]]))
    
    ## check whether conditional permutation is useful
    if(nlevels(CondPartitions) == nlevels(AllPartitions)) return(NULL)
    
    ## select only the unique partitions. Permuting the partitions
    ## that are the same after including the variable of interest 
    ## cannot change the prediction, and therefore is redundant.
    PartsWithout <- listParts(CondPartitions, oob)
    PartsWith <- listParts(AllPartitions, oob)
    parts <- GetUniqueParts(PartsWithout, PartsWith)
  }
  return(vector_lengths(parts))
} 





vector_lengths <- function(vector_list) {
  # Initialize an empty vector to store the lengths
  lengths <- numeric(length(vector_list))
  
  # Iterate through the list and calculate the length of each vector
  for (i in 1:length(vector_list)) {
    lengths[i] <- length(vector_list[[i]])
  }
  
  # Return the vector of lengths
  return(lengths)
}
