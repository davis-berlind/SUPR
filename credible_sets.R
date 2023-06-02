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

cred_set_prisca <- function(probs, max_length = floor(0.5 * nrow(probs)), level = 0.9, merge_prob = 0.5) {
  L <- ncol(probs)
  T <- nrow(probs)
  cs <- list()
  skip <- c()
  for(l in 1:L) {
    if (l %in% skip) next
    set <- order(probs[,l], decreasing = TRUE)[1:which.max(cumsum(probs[order(probs[,l], decreasing = TRUE), l]) > level)]
    if (l < L) {
      for (k in (l+1):L) {
        if (sum(probs[,l] * probs[,k] >= merge_prob)) {
          set <- union(set, order(probs[,k], decreasing = TRUE)[1:which.max(cumsum(probs[order(probs[,k], decreasing = TRUE), k]) > level)])
          skip <- c(skip, k)
        }
      }
    }
    if (length(set) <= max_length) {
      cs <- c(cs,list(set[order(set)]))
    }
  }
  
  if (length(cs) == 0) return(cs)
  
  cs <- unique(cs)
  # nset <- length(cs)
  # subset <- c()
  # for (i in 1:nset) {
  #   for(j in 1:nset) {
  #     if (i == j) next
  #     if (length(setdiff(cs[[i]], cs[[j]])) == 0) subset <- c(subset, i)
  #   }
  # }
  # if (length(subset) > 0) cs[[subset]] <- NULL
  
  return(cs)
}
