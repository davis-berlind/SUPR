library(latex2exp)
set.seed(3)

T <- 1000
J <- 10
C <- 100
min_space <- 50
cp <- hsmuce_simulation(T, J, C, min_space)

fit <- mich(cp$y, J = 10, verbose = TRUE, tol = 1e-7, 
            fit_intercept = FALSE, fit_scale = FALSE)
est_cp <- mich_sets(fit$meanvar_model$pi_bar)$cp

fit_rev <- mich(cp$y, J = 10, verbose = TRUE, tol = 1e-7, reverse = TRUE)
est_cp_rev <- mich_sets(fit_rev$meanvar_model$pi_bar)$cp

par(mfrow = c(2,1), oma = c(0,0,0,0), mar = c(4,4,1,2))

plot(cp$y, type = "l", ylab = "", xlab = "", xaxt = "n")
axis(side = 1, at = cp$changepoints, lwd.ticks = 2, labels = FALSE)
lines(fit$mu, col = "red", lwd = 3)
lines(fit$mu + 2 / sqrt(fit$lambda), col = "red", lwd = 1, lty = 2)
lines(fit$mu - 2 / sqrt(fit$lambda), col = "red", lwd = 1, lty = 2)
abline(v = est_cp, col = "red", lwd = 1, lty = 2)

plot(cp$y, type = "l", ylab = "", xlab = "Index", xaxt = "n")
axis(side = 1, at = cp$changepoints, lwd.ticks = 2, labels = TRUE)
lines(fit_rev$mu, col = "blue", lwd = 3)
lines(fit_rev$mu + 2 / sqrt(fit_rev$lambda), col = "blue", lwd = 1, lty = 2)
lines(fit_rev$mu - 2 / sqrt(fit_rev$lambda), col = "blue", lwd = 1, lty = 2)
abline(v = est_cp_rev, col = "blue", lwd = 1, lty = 2)

