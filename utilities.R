#### helper functions #### 

revcumsum <- function(x){
  rev(cumsum(rev(x)))
}

lambda_bar_fn <- function(u, v, prob) {
  cumsum((u / v) * prob) + c(revcumsum(prob[-1]), 0)
}

#### data formatting functions #### 

prob_check <- function(probs, fit.intercept, n, T) {
  test <- TRUE
  if (is.null(probs)) {
    probs <- matrix(1/T, ncol = n, nrow = T)
  } else {
    # delete prob that first point is a change if fitting intercept and too many pi
    if (fit.intercept & nrow(probs) == (T + 1)) probs <- probs[-1,] / colSums(probs[-1,])  
    if (!is.numeric(probs)) test <- FALSE
    if (is.array(probs) & (nrow(probs) != T | ncol(probs) != n | any(round(colSums(probs), 10) != 1))) test <- FALSE
  }
  if (test) return(probs)
  else return(NULL)
}

prior_check <- function(prior, n) {
  test <- TRUE
  if (!is.numeric(prior)) test <- FALSE
  else if (any(prior < 0) | (length(prior) > 1 & length(prior) != n)) test <- FALSE
  if (!test) return(NULL)
  else if (length(prior) == 1) return(rep(prior, n))
  else return (prior)
}
