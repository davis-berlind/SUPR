
# Data Generation

set.seed(10)

T <- 800

beta <- rep(0, T)
beta[201:400] <- 3
beta[401:600] <- -2
beta[601:T] <- -4

lambda <- rep(3, T)
lambda[201:400] <- 10
lambda[401:600] <- .5
lambda[601:T] <- 2

y <- rnorm(T, mean = beta, sd = sqrt(1 / lambda))

diags <- c()
for(i in 1:10) {
  mich_fit <- mich_ii(y, L = i, fit.intercept = TRUE)
  diags <- c(diags, sum(diag(round(t(mich_fit$pi) %*% mich_fit$pi,2))))
}

plot_mich(mich_fit)

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
