
# Data Generation

set.seed(10)

T <- 600
X <- matrix(0, nrow = T, ncol = T)
X[lower.tri(X, diag = TRUE)] <- 1

beta <- rep(0, T)
beta[201] <- 3
beta[401] <- -5

lambda <- rep(1, T)
lambda[301:T] <- 10

y <- rnorm(T, mean = X %*% beta, sd = sqrt(1 / lambda))

supr_mod <- supr(y, sigma2 = 1, K = 5, L = 5, sigma2_0 = 0.1)


plot(y)
lines(rowSums(X %*% supr_mod$beta), col = "red", lwd = 3)
lines(rowSums(X %*% supr_mod$beta) + 2 / supr_mod$lambda2, lty = 3, col = "blue", lwd = 3)
lines(rowSums(X %*% supr_mod$beta) - 2 / supr_mod$lambda2, lty = 3, col = "blue", lwd = 3)

# pips
plot(1 - apply(1 - supr_mod$gamma, 1, prod))
abline(v = 201)
abline(v = 401)

plot(1 - apply(1 - supr_mod$alpha, 1, prod))
abline(v = 201)
abline(v = 401)
