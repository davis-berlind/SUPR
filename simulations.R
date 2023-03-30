set.seed(7)

# simultion parameters
T <- c(200, 500, 1000)
N <- 100
loc_sd <- sqrt(10) # same order as (exp((log(10) / 2)^2) - 1) * exp((log(10) / 2)^2)
scale_sd <- log(10) / 2 # justify this 

# supr parameters
u_0 <- 0.001
v_0 <- 0.001
sigma2_0 <- 10  # large enough?
tol <- 1e-5 # max iter + relax tol 



results <- data.frame(matrix(0,nrow = 0, ncol = 13))

bias_loc <- rep(0, N) 
bias_scl <- rep(0, N) 
haus_loc <- rep(0, N) 
haus_scl <- rep(0, N) 
length_loc <- rep(0, N) 
length_scl <- rep(0, N) 
precision_loc <- rep(0, N) 
precision_scl <- rep(0, N) 
recall_loc <- rep(0, N) 
recall_scl <- rep(0, N) 
f_loc <- rep(0, N) 
f_scl <- rep(0, N) 

for (t in T) {
  L_true <- floor(sqrt(t) / 4)
  K_true <- L_true
  
  L <- floor(t / 30)
  K <- L
  
  for (i in 1:N) {
    active_loc <- point_picker(t, L_true)
    active_scl <- point_picker(t, K_true)
    
    mn <- c(0, rnorm(L_true, sd = loc_sd))
    mn <- cumsum(mn)
    mn <- mn[as.numeric(cut(1:t, breaks = c(0, active_loc, t)))]    
    
    scl <- c(1, exp(rnorm(K_true, sd = scale_sd)))
    scl <- scl[as.numeric(cut(1:t, breaks = c(0, active_scl, t)))]    
    
    y <- rnorm(t, mn, sqrt(scl))
    
    model <- supr(y, sigma2 = 1, sigma2_0, K, L, tol, u_0 , v_0)
    
    loc_cs <- cred_set_susie(model$gamma)
    scl_cs <- cred_set_prisca(model$alpha)
    
    bias_loc[i] <- L_true - length(loc_cs)
    bias_scl[i] <- K_true - length(scl_cs)
    
    loc_cs <- ifelse(length(loc_cs) == 0, 0, loc_cs)
    scl_cs <- ifelse(length(scl_cs) == 0, 0, scl_cs)
    
    haus_loc[i] <- ifelse(length(unlist(loc_cs)) == 1,
                          min(abs(sapply(active_loc, rep, length = length(unlist(loc_cs))) - unlist(loc_cs))),
                          max(apply(abs(sapply(active_loc, rep, length = length(unlist(loc_cs))) - unlist(loc_cs)), 2, min)))
    haus_scl[i] <- ifelse(length(unlist(scl_cs)) == 1,
                          min(abs(sapply(active_scl, rep, length = length(unlist(scl_cs))) - unlist(scl_cs))),
                          max(apply(abs(sapply(active_scl, rep, length = length(unlist(scl_cs))) - unlist(scl_cs)), 2, min)))
    length_loc[i] <- mean(sapply(loc_cs, length))
    length_scl[i] <- mean(sapply(scl_cs, length))
    precision_loc[i] <- sum(active_loc %in% unlist(loc_cs)) / length(loc_cs)
    precision_scl[i] <-  sum(active_scl %in% unlist(scl_cs)) / length(scl_cs)
    recall_loc[i] <- mean(active_loc %in% unlist(loc_cs))
    recall_scl[i] <- mean(active_scl %in% unlist(scl_cs))
    f_loc[i] <- 2 * precision_loc[i] * recall_loc[i] / (precision_loc[i] + recall_loc[i])
    f_scl[i] <- 2 * precision_scl[i] * recall_scl[i] / (precision_scl[i] + recall_scl[i])
    
    print(i)
  }
  results <- rbind(results, 
                   c(t, mean(bias_loc), mean(bias_scl),
                     mean(haus_loc), mean(haus_scl),
                     mean(length_loc), mean(length_scl), 
                     mean(precision_loc), mean(precision_scl),
                     mean(recall_loc), mean(recall_scl), 
                     mean(f_loc), mean(f_scl)))
}

names(results) <- c("T", "bias_loc", "bias_scl", "haus_loc", "haus_scl",
                    "length_loc", "length_scl", "precision_loc", "precision_scl",
                    "recall_loc", "recall_scl", "f_loc", "f_scl")

point_picker <- function(T, K) {
  space <- min(floor(sqrt(T)), 30)
  valid <- 2:(T-2)
  picked <- c()
  for (k in 1:K) {
    pick <- sample(valid, size = 1)
    picked <- c(picked, pick)
    valid <- valid[valid < pick - space | valid > pick + space]
  }
  return(picked)
}
