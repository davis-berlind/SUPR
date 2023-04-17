
# Data Generation

set.seed(10)

T <- 800

beta <- rep(1, T)
beta[201:400] <- 3
beta[401:600] <- -2
beta[601:T] <- -4

lambda <- rep(2, T)
lambda[201:400] <- 10
lambda[401:T] <- .5
#lambda[601:T] <- 2

y <- rnorm(T, mean = beta, sd = sqrt(1 / lambda))

supr_fit <- supr(y, K = 5, L = 5)

plot_supr_fit(supr_fit)

cred_set_prisca(supr_fit$pi)
cred_set_prisca(supr_fit$omega)
cred_set_susie(supr_fit$pi)
cred_set_susie(supr_fit$omega)



supr_mod <- supr_sim(y, sigma2 = 1, sigma2_0 = 10, L = 5, tol = 1e-10, u_0 = 1e-3, v_0 = 1e-3)

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

cred_set_susie(supr_mod$pi)
