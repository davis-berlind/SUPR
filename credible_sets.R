#### credible set functions  #### 

cred_set_susie <- function(probs, fit.intercept = TRUE, level = 0.95) {
  L <- ncol(probs)
  T <- nrow(probs) + fit.intercept
  cs <- list()
  for(l in 1:L) {
    set <- order(probs[,l], decreasing = TRUE)[1:which.max(cumsum(probs[order(probs[,l], decreasing = TRUE), l]) > level)]
    set <- set[order(set)]
    if(length(set) == 1) {
      cs <- c(cs, list(set))
      next
    } 
    cmbs <- combn(set, 2)
    i <- cmbs[1,] + fit.intercept
    j <- cmbs[2,] + fit.intercept
    E_i <- 1 - (i - 1) / T  
    E_j <- 1 - (j - 1) / T
    purity <- min(sqrt(E_j) * (1 - E_i) / sqrt(E_i * (1 - E_i) * (1 - E_j)))
    if (purity > 0.5) {
      cs <- c(cs,list(set[order(set)]))
    }
  }
  
  if (length(cs) == 0) return(cs)
  
  cs <- unique(cs)
  nset <- length(cs)
  subset <- c()
  for (i in 1:nset) {
    for(j in 1:nset) {
      if (i == j) next
      if (length(setdiff(cs[[i]], cs[[j]])) == 0) subset <- c(subset, i)
    }
  }
  
  if (length(subset) > 0) cs[[subset]] <- NULL
  
  return(cs)
}

#' Posterior Credible Set
#' 
#' Given posterior probabilities for the location of a change-point, this 
#' function finds the fewest possible locations needed to construct a set with
#' a posterior probability of containing the change-point that greater than or
#' equal to some specified threshold.
#'
#' @param prob  A numeric vector. A vector of posterior probabilities for the 
#'              location of the change-point.
#' @param level A scalar. A single number in (0,1) that gives the lower bound
#'              for the probability that the credible set contains a 
#'              change-point 
#'
#' @return Level `level` posterior credible set for the location of a 
#'         change-point. 
#'         
cred_set <- function(prob, level) {
  order(prob, decreasing = TRUE)[1:which.max(cumsum(prob[order(prob, decreasing = TRUE)]) > level)]
}

mich_sets <- function(probs, max_length = log(nrow(probs))^2 * 2, level = 0.9) {
  
  # initialize credible sets and changepoints
  cs <- list()
  est_cp <- numeric(0)
  
  cred_sets <- apply(probs, 2, cred_set, level = level, simplify = FALSE)
  
  # drop probs/sets that are longer than max_length
  keep <- which(sapply(cred_sets, length) <= max_length)
  if (length(keep) > 0) {
    est_cp <- apply(probs[, keep, drop = FALSE], 2, which.max)
    cs <- lapply(keep, function(i) cred_sets[[i]])
    # order change-points
    cs <- lapply(order(est_cp), function(i) cs[[i]])
    est_cp <- est_cp[order(est_cp)]
  }

  return(list(cp = est_cp, sets = sapply(cs, function(set) set[order(set)])))
}

