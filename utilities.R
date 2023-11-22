
revcumsum <- function(x){
  T <- length(x)
  s <- c(0, cumsum(x[-T]))
  return(x[T] + s[T] - s)
}

prob_check <- function(probs, n, T) {
  if (is.null(probs)) return(matrix(nrow = T, ncol = n))
  if (is.vector(probs)) probs <- sapply(1:max(1,n), function(i) probs)
  if (is.array(probs)) {
    if (nrow(probs) != T) return(NULL)
    if (any(round(colSums(probs), 10) != 1)) return(NULL)
    if (n > 0 & ncol(probs) != n) return(NULL)
  }
  return(probs)
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

integer_check <- function(n) {
  if (is.numeric(n)) {
    if  (n == round(n) & n >= 0) return (n)
  } else return(NULL)
}
