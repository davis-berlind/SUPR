library(latex2exp)
set.seed(1)

T <- 1000
J <- 10
C <- 100
min_space <- 50
cp <- hsmuce_simulation(T, J, C, min_space)

fit <- mich(cp$y, J = 10, merge_prob = 1, verbose = TRUE, tol = 1e-5, 
            fit_intercept = FALSE, fit_scale = FALSE)
est_cp <- mich_sets(fit$meanvar_model$pi_bar)$cp

fit_merge <- mich(cp$y, J = 10, verbose = TRUE, tol = 1e-5, 
                  fit_intercept = FALSE, fit_scale = FALSE)
est_cp_merge <- mich_sets(fit_merge$meanvar_model$pi_bar)$cp

# merge probs
round(t(fit$meanvar_model$pi_bar) %*% fit$meanvar_model$pi_bar,3)

par(mfrow = c(2,1), oma = c(0,0,0,0), mar = c(4,4,1,2))

plot(cp$y, type = "l", ylab = "y", xlab = "", xaxt = "n")
axis(side = 1, at = cp$changepoints, lwd.ticks = 2)
lines(fit$mu, col = "red", lwd = 3)
lines(fit$mu + 2 / sqrt(fit$lambda), col = "red", lwd = 1, lty = 2)
lines(fit$mu - 2 / sqrt(fit$lambda), col = "red", lwd = 1, lty = 2)
abline(v = est_cp, col = "red", lwd = 1, lty = 2)

# legend("topleft", legend = c(TeX("\\beta = 1"),  TeX("$\\beta = \\log^2 T / T^2$")), 
#        col = c("red", "blue"), lty = c(1,2), lwd = 3, cex = 1.25)

dex = 170:210
barplot(t(fit$meanvar_model$pi_bar[dex,c(1,9)]), beside = TRUE, col = c("red", "darkred"), 
        names.arg = dex, cex.names = 0.75, xlab = "Index", ylab = "Posterior Probability")
legend("topleft", legend = c(TeX("$Pr(\\tau_2 = t\\;|\\;y_{1;T})$"),  TeX("$Pr(\\tau_3 = t\\;|\\;y_{1;T})$")), 
       fill = c("red", "darkred"), cex = 1)

dev.off()
par(oma = c(0,0,0,0), mar = c(4,4,1,2))
plot(cp$y, type = "l", ylab = "y", xlab = "Index", xaxt = "n")
axis(side = 1, at = cp$changepoints, lwd.ticks = 2)
lines(fit_merge$mu, col = "blue", lwd = 3)
lines(fit_merge$mu + 2 / sqrt(fit_merge$lambda), col = "blue", lwd = 1, lty = 2)
lines(fit_merge$mu - 2 / sqrt(fit_merge$lambda), col = "blue", lwd = 1, lty = 2)
abline(v = est_cp_merge, col = "blue", lwd = 1, lty = 2)

