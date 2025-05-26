set.seed(222)

# mean-scp ####
y <- c(rnorm(50), rnorm(50, 1))

fit <- mich(y, L = 1, fit_intercept = FALSE, fit_scale = FALSE)
cp <- mich_sets(fit$mean_model$pi_bar)$cp
cs <- mich_sets(fit$mean_model$pi_bar)$sets

png("~/Desktop/mean_scp.png", width = 1350, height = 850)

plot(y, type = "l", lwd = 2, main = "Mean-SCP Model",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(v = 50, col = "red", lty = 2, lwd = 2)
dev.off()

png("~/Desktop/mean_scp_post.png", width = 1350, height = 850)

plot(y, type = "l", lwd = 2, main = "Mean-SCP Fit with 90% Credible Set",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)

left = min(cs[[1]])
key <- diff(cs[[1]])
for(i in 1:(length(cs[[1]]) - 1)) {
  if (key[i] != 1) {
    rect(xleft = left-0.5, xright = cs[[1]][i] + 0.5, 
         ybottom = par("usr")[3], ytop = par("usr")[4], 
         col =  adjustcolor("lightblue", alpha = 0.3), 
         border = NA)
    left = cs[[1]][i+1]
  }
}
abline(v = cp, col = "blue", lty = 2, lwd = 2)
abline(v = 50, col = "red", lty = 2, lwd = 2)

# lines(c(rep(0,50), rep(1,50)), col = "red", lwd = 2)
# lines(fit$mu, col = "blue", lwd = 2)
legend("bottomright", inset = 0.01, legend = c("True CP", "Estimated CP"), 
       col = c("red", "blue") , lty = 2, lwd = 2, cex = 1.5)

dev.off()

# var-scp ####
y <- c(rnorm(50), rnorm(50, 0, 0.5))

fit <- mich(y, K = 1, fit_intercept = FALSE, fit_scale = FALSE)
cp <- mich_sets(fit$var_model$pi_bar)$cp
cs <- mich_sets(fit$var_model$pi_bar)$sets

png("~/Desktop/var_scp.png", width = 1350, height = 850)

plot(y, type = "l", lwd = 2, main = "Var-SCP Model",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(v = 50, col = "red", lty = 2, lwd = 2)

dev.off()

png("~/Desktop/var_scp_post.png", width = 1350, height = 850)

plot(y, type = "l", lwd = 2, main = "Var-SCP Fit with 90% Credible Set",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)

left = min(cs[[1]])
key <- diff(cs[[1]])
for(i in 1:(length(cs[[1]]) - 1)) {
  if (key[i] != 1) {
    rect(xleft = left-0.5, xright = cs[[1]][i] + 0.5, 
         ybottom = par("usr")[3], ytop = par("usr")[4], 
         col =  adjustcolor("lightblue", alpha = 0.3), 
         border = NA)
    left = cs[[1]][i+1]
  }
}
abline(v = cp, col = "blue", lty = 1, lwd = 2)
abline(v = 50, col = "red", lty = 2, lwd = 2)

# lines(1.96 * 1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 1)
# lines(-1.96 * 1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 1)
# 
# lines(c(rep(1.96,50), rep(1.95*0.5,50)), col = "red", lwd = 2, lty = 2)
# lines(-c(rep(1.96,50), rep(1.95*0.5,50)), col = "red", lwd = 2, lty = 2)
legend("bottomright", inset = 0.01, legend = c("True CP", "Estimated CP"), 
       col = c("red", "blue") , lty = 2, lwd = 2, cex = 1.5)

dev.off()

# meanvar-scp ####
y <- c(rnorm(50), rnorm(50, 0.5, 0.5))

fit <- mich(y, J = 1, fit_intercept = FALSE, fit_scale = FALSE)
cp <- mich_sets(fit$meanvar_model$pi_bar)$cp
cs <- mich_sets(fit$meanvar_model$pi_bar)$sets

png("~/Desktop/meanvar_scp.png", width = 1350, height = 850)

plot(y, type = "l", lwd = 2, main = "MeanVar-SCP Model",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
abline(v = 50, col = "red", lty = 2, lwd = 2)

dev.off()

png("~/Desktop/meanvar_scp_post.png", width = 1350, height = 850)

plot(y, type = "l", lwd = 2, main = "MeanVar-SCP Fit with 90% Credible Set",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)

left = min(cs[[1]])
key <- diff(cs[[1]])
for(i in 1:(length(cs[[1]]) - 1)) {
  if (key[i] != 1) {
    rect(xleft = left-0.5, xright = cs[[1]][i] + 0.5, 
         ybottom = par("usr")[3], ytop = par("usr")[4], 
         col =  adjustcolor("lightblue", alpha = 0.3), 
         border = NA)
    left = cs[[1]][i+1]
  }
}
abline(v = cp, col = "blue", lty = 2, lwd = 2)
abline(v = 50, col = "red", lty = 2, lwd = 2)

# lines(fit$mu, col = "blue", lwd = 2, lty = 1)
# lines(fit$mu + 1.96 * 1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 1)
# lines(fit$mu -1.96 * 1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 1)
# 
# lines(c(rep(0,50), rep(0.5,50)), col = "red", lwd = 2, lty = 2)
# lines(c(rep(1.96,50), 0.5 + rep(1.95*0.5,50)), col = "red", lwd = 2, lty = 2)
# lines(c(rep(-1.96,50),0.5 - rep(1.95*0.5,50)), col = "red", lwd = 2, lty = 2)
legend("bottomright", inset = 0.01, legend = c("True CP", "Estimated CP"), 
       col = c("red", "blue") , lty = 2, lwd = 2, cex = 1.5)

dev.off()