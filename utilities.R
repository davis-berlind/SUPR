
revcumsum <- function(x){
  T <- length(x)
  s <- c(0, cumsum(x[-T]))
  return(x[T] + s[T] - s)
}

prob_check <- function(probs, n, T) {
  if (length(probs) == 1) {
    if (probs %in% c("uniform", "weighted")) return(probs)
  } else if (is.vector(probs) & is.numeric(probs)) {
    if (all(!is.na(probs)) & (length(probs) == T) & is.numeric(probs)) {
      if (round(sum(probs), 10) == 1) return(sapply(1:max(1, n), function(i) probs))
    }
  } else if (is.array(probs) & is.numeric(probs)) {
    if (all(!is.na(probs)) & (nrow(probs) == T & ncol(probs) == n)) {
      if (all(round(colSums(probs), 10) == 1))  return(probs)
    }
  }
  stop(paste0(deparse(substitute(probs)), 
              " must be 'uniform', 'weighted', or a length T vector or a T x ",
              deparse(substitute(n)),
              " matrix with columns that sum to one."))
}

logical_check <- function(x) {
  if (is.logical(x) & length(x) == 1) {
    if (!is.na(x)) return(x) 
  }
  stop(paste0(deparse(substitute(x)), " must be TRUE or FALSE."))
}

scalar_check <- function(x) {
  if (is.numeric(x) & length(x) == 1) {
    if (!is.na(x)) {
      if (x > 0) return(x) 
    }
  } 
  stop(paste0(deparse(substitute(x)), " must be an scalar > 0."))
}

integer_check <- function(n) {
  if (is.numeric(n) & length(n) == 1) {
    if (!is.na(n)) {
      if (n == round(n) & n >= 0) return(n) 
    }
  }
  stop(paste0(deparse(substitute(n)), " must be an integer >= 0."))
}
