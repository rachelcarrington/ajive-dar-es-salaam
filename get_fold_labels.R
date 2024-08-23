### Split data into k folds for cross-validation

get_fold_labels <- function(k, n){

  ### Divide n data points randomly into k even folds
  
  fold_size <- floor(n / k)
  remainder <- n - fold_size * k
  fold_labels <- c(rep(1:k, fold_size), 1:remainder)
  fold_labels <- sample(fold_labels, n) ## randomly re-order labels

  return(fold_labels)
}
