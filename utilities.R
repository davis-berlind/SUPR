#### helper functions #### 

revcumsum <- function(x){
  T <- length(x)
  s <- c(0, cumsum(x[-T]))
  return(x[T] + s[T] - s)
}

#### data formatting functions #### 

prob_check <- function(probs, n, T) {
  if (n == 0) return (matrix(ncol = n, nrow = T))
  test <- TRUE
  if (is.null(probs)) {
    probs <- matrix(1 / T, ncol = n, nrow = T)
  } else {
    if (!is.numeric(probs)) test <- FALSE
    if (is.array(probs) & (nrow(probs) != T | ncol(probs) != n | any(round(colSums(probs), 10) != 1))) test <- FALSE
  }
  if (test) return(probs)
  else return(NULL)
}

prior_check <- function(prior, n) {
  if (length(n) == 0) return(numeric())
  test <- TRUE
  if (!is.numeric(prior)) test <- FALSE
  else if (any(prior < 0) | (length(prior) > 1 & length(prior) != n)) test <- FALSE
  if (!test) return(NULL)
  else if (length(prior) == 1) return(rep(prior, n))
  else return (prior)
}

component_check <- function(n) {
  if (is.character(n)) {
    if (n == "auto") return(n)
  } else if (is.numeric(n)) {
    if  (n == round(n) & n >= 0) return (n)
  } else return(NULL)
}
