L = 3; J = 0; K=3; fit.intercept = TRUE; fit.scale = TRUE;
tol = 1e-5; B_l = 1; B_r = B_l; verbose = FALSE;
tau_j = 0.1; tau_l = tau_j; 
u_j = 1e-3; v_j = u_j; u_k = u_j; v_k = v_j;
theta = NULL; pi = NULL; omega = NULL;
# Data Generation

set.seed(10)

T <- 800

mu <- rep(0, T)
mu[201:400] <- 3
mu[401:600] <- -2
mu[601:T] <- -4

lambda <- rep(1, T)
lambda[201:400] <- 10
lambda[401:600] <- .5
lambda[601:T] <- 2

y <- rnorm(T, mean = mu, sd = sqrt(1 / lambda))
tmp=mich(y, J = 3, verbose = TRUE, fit.intercept = FALSE, fit.scale = FALSE)
tmp=mich(y, K = 3, L = 3, verbose = TRUE, fit.intercept = TRUE, fit.scale = TRUE)

tmp=mich(dat$y, J=4, L=4, K=4, verbose = TRUE, fit.intercept = FALSE, fit.scale = FALSE, conv_crit = "ELBO", tol = 1e-3)
plot(dat$y)
abline(v=apply(tmp$mean.scale.model$probs,2,which.max),col="blue")
abline(v=apply(tmp$mean.model$probs,2,which.max),col="red",lty=2)
abline(v=apply(tmp$scale.model$probs,2,which.max),col="orange",lty=3)
lines(tmp$mu)
plot(tmp$elbo)

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
