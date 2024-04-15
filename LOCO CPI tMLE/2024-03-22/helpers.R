### Helper Functions
first_term_noresidue_rf <-
  function(index, full_data, full_model, reduced_model) {
    mean(2*(full_model$predictions[index] 
            - reduced_model$predictions[index]))
  }

first_term_noresidue_cv <-
  function(index, full_data, full_model, reduced_model) {
    mean(2*(full_model$predictions[index] 
            - reduced_model$predictions[index]))
  }

first_term_noresidue_gam <- 
  function(index, full_data, full_model, reduced_model) {
    mean(2*(full_model$fitted.values[index] 
            - reduced_model$fitted.values[index]))
  }

first_term <-
  function(index, full_data, full_model) {
    mean(2*(full_data[index, 1] 
            - full_model$predictions[index])*(full_model$predictions[index] 
                                              - reduced_model$predictions[index]))
  }


logit <- function(p) {
  log(p / (1 - p))
}

rescaling <- function(curr_frame, covariate_index){
  dist = max(curr_frame[,1]) - min(curr_frame[,1])
  min  = min(curr_frame[,1])
  Y_trunc = (curr_frame[,1] - min)/dist
  output.fullframe = data.frame(y = Y_trunc, X = curr_frame[,-1])
  colnames(output.fullframe) <- c("output", 1:10)
  output.reducedframe = data.frame(y = Y_trunc, X = curr_frame[,-c(1, covariate_index+1)])
  colnames(output.reducedframe) <- c("output", 2:10)
  return(list(output.fullframe, output.reducedframe, Y_trunc, c(dist, min)))
}

rescaling_cv <- function(curr_frame, covariate_index, dist, min){
  Y_trunc = (curr_frame[,1] - min)/dist
  output.fullframe = data.frame(y = Y_trunc, X = curr_frame[,-1])
  colnames(output.fullframe) <- c("output", 1:10)
  output.reducedframe = data.frame(y = Y_trunc, X = curr_frame[,-c(1, covariate_index+1)])
  colnames(output.reducedframe) <- c("output", 2:10)
  return(list(output.fullframe, output.reducedframe, Y_trunc, c(dist, min)))
}


trunc <- function(full_model_pred, reduced_model_pred){
  trunc_index = full_model_pred>.005&full_model_pred<(1-.005)
  trunc_reduce =  reduced_model_pred>.005&reduced_model_pred<(1-.005)
  trunc_index = as.logical(trunc_index*trunc_reduce)
  return(trunc_index)
}










