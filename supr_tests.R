
# Data Generation

set.seed(10)

T <- 800
X <- matrix(0, nrow = T, ncol = T)
X[lower.tri(X, diag = TRUE)] <- 1

beta <- rep(0, T)
beta[201] <- 3
beta[401] <- -2
beta[600] <- -4

lambda <- rep(1, T)
lambda[201:400] <- 10
lambda[401:T] <- .5
lambda[601:T] <- 2

y <- rnorm(T, mean = X %*% beta, sd = sqrt(1 / lambda))

supr_mod <- supr(y, sigma2 = 1, K = 3, L = 3, sigma2_0 = 10)
supr_mod <- supr_sim(y, sigma2 = 1, sigma2_0 = 10, L = 3, tol = 1e-5, u_0 = 1e-3, v_0 = 1e-3)

plot(y)
lines(rowSums(X %*% supr_mod$beta), col = "red", lwd = 3)
lines(rowSums(X %*% supr_mod$beta) + 2 / sqrt(supr_mod$lambda2), lty = 3, col = "blue", lwd = 3)
lines(rowSums(X %*% supr_mod$beta) - 2 / sqrt(supr_mod$lambda2), lty = 3, col = "blue", lwd = 3)

# pips
plot(1 - apply(1 - supr_mod$gamma, 1, prod))
abline(v = 201)
abline(v = 401)

plot(1 - apply(1 - supr_mod$alpha, 1, prod))
abline(v = 201)
abline(v = 401)
