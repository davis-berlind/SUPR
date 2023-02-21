set.seed(10)

# Data Generation
T <- 100
X <- matrix(0, nrow = T, ncol = T)
X[lower.tri(X, diag = TRUE)] <- 1

beta <- rep(0, T)
beta[34] <- 3
beta[67] <- -5

lambda <- rep(1, T)
lambda[25:50] <- 10
#lambda[51:T] <- 0.33

y <- rnorm(T, mean = X %*% beta, sd = sqrt(1 / lambda))

# function
tol <- 1e-7
T <- length(y)
X <- matrix(0, nrow = T, ncol = T)
X[lower.tri(X, diag = TRUE)] <- 1
  
sigma2 <- 1
sigma2_0 <- 1
u_0 <- 1e-3
v_0 <- 1e-3
pi <- 1 / T
omega <- 1 / T

L <- 5
K <- 5

# initializing parameters

gamma_post <- matrix(1 / T, nrow = T, ncol = L)
alpha_post <- matrix(1 / T, nrow = T, ncol = K)

mu_post <- matrix(0, nrow = T, ncol = L)
sigma_post <- rep(1, length = T)
beta <- mu_post * gamma_post
u_post <- u_0 + (T - 1:T + 1) / 2
v_post <- matrix(u_post, nrow = T, ncol = K)

lambda2_k <- matrix(0, nrow = T, ncol = K)

for (k in 1:K) {
  tmp <- t(u_post / v_post[,k] * t(X))
  tmp[upper.tri(tmp)] <- 1
  lambda2_k[,k] <- tmp %*% alpha_post[,k]
}

lambda2 <- apply(lambda2_k, 1, prod)

params <- c(sigma_post, mu_post, gamma_post, v_post, alpha_post)

while (TRUE) {
  
  new_params <- c()
  
  r_bar <- y - rowSums(X %*% beta)
  sigma_post <- 1 / (revcumsum(lambda2) / sigma2 + 1 / sigma2_0)
  
  new_params <- c(new_params, sigma_post)
  
  for (l in 1:L) {

    r_bar <- r_bar + X %*% beta[,l]
    mu_post[,l] <- (sigma_post / sigma2) * revcumsum(lambda2 * r_bar)

    gamma_post[,l] <- log(pi) + 0.5 * log(sigma_post) + mu_post[,l]^2 / (2 * sigma_post)
    gamma_post[,l] <- prop.table(exp(gamma_post[,l] - max(gamma_post[,l])))

    beta[,l] <- mu_post[,l] * gamma_post[,l]

    r_bar <- r_bar - X %*% beta[,l]
  }
  
  new_params <- c(new_params, mu_post, gamma_post)
  
  for (k in 1:K) {
    lambda2 <- lambda2 / lambda2_k[,k]
    
    z2_tilde <- lambda2 * ((y - rowSums(X %*% beta))^2 - rowSums((X %*% beta)^2) + cumsum(rowSums((mu_post^2 + sigma_post) * gamma_post)))
    
    v_post[,k] <- v_0 + revcumsum(z2_tilde) / (2 * sigma2)
    
    alpha_post[,k] <- log(omega) + lgamma(u_post) - u_post * log(v_post[,k]) - cumsum(c(0,z2_tilde))[-(T+1)] / (2 * sigma2)
    alpha_post[,k] <- prop.table(exp(alpha_post[,k] - max(alpha_post[,k])))
    
    tmp <- t(u_post / v_post[,k] * t(X))
    tmp[upper.tri(tmp)] <- 1
    lambda2_k[,k] <- tmp %*% alpha_post[,k]
    
    lambda2 <- lambda2 * lambda2_k[,k]
  }
  
  new_params <- c(new_params, v_post, alpha_post)
  
  if (sum((params - new_params)^2) < tol) break
  params <- new_params
}

revcumsum <- function(x){
  rev(cumsum(rev(x)))
}

plot(y, ylim = c(min(rowSums(X%*%beta) - 2 * 1 / lambda2), max(rowSums(X%*%beta) + 2 * 1 / lambda2)))
lines(rowSums(X%*%beta), col = "red")
lines(rowSums(X%*%beta) + 2/lambda2, lty = 3)
lines(rowSums(X%*%beta) - 2/lambda2, lty = 3)

# pips
plot(1 - apply(1 - gamma_post, 1, prod))
abline(v = 34)
abline(v = 67)

plot(1 - apply(1 - alpha_post, 1, prod))
abline(v = 51)
abline(v = 25)

cred_set <- function(probs, level = 0.9) {
  L <- ncol(probs) 
  cs <- list()
  for(l in 1:L) {
    set <- order(probs[,l],decreasing = TRUE)[1:which.max(cumsum(probs[order(probs[,l], decreasing = TRUE), l]) > level)]
    if(length(set) == 1) {
      cs[[l]] <- set
      next
    }
    combs <- combn(set,2)
    p_i <- probs[combs[1,],l]
    p_j <- probs[combs[2,],l]
    if(min(abs(p_i * p_j / (sqrt(p_i*(1-p_i) *p_j*(1-p_j))))) > 0.5) {
      cs[[l]] <- set
    }
  }
  return(cs)
}

cred_set(alpha_post)
cred_set(gamma_post)
