ion <- read.csv("~/Desktop/ICdata.csv")
id <- seq(403300, 435811, 11)
x <- ion$time[id]
unit <- x[2] - x[1]
y <- ion$data[id]
T <- length(y)

dev.off()
par(oma = c(0,1,0,0))

plot(x, y, type = "l", main = "", 
     xlab = "Time (s)", ylab = "Conductance (nS)", 
     cex.main=2, cex.lab=1.5)

library(stepR)
library(changepoint)
library(latex2exp)

png("~/Desktop/ion_plot_1.png", width = 1250, height = (1250 / 2) * 3/2)
par(mfrow = c(3,1), oma = c(0,1,0,0), mar = c(4,4,3,2))

#### MICH MeanVar ####
fit_mean_var <- mich(y, J_auto = TRUE, max_iter = Inf, verbose = TRUE, tol = 1e-6)
cred_sets <- mich_sets(fit_mean_var$meanvar_model$pi_bar, level = 0.95)
mu <- fit_mean_var$mu
est_cp <- cred_sets$cp
plot(x, y, type = "l", main = paste0("MICH (J = ", fit_mean_var$J,")"), 
     xlab = "", ylab = "", col = "white", 
     cex.main=2.5, cex.lab=2, cex.axis=1.75, )
for(i in unlist(cred_sets$sets)) { 
  rect(xleft = x[i]-2*unit, xright = x[i] + 2*unit, 
       ybottom = par("usr")[3], ytop = par("usr")[4], 
       col =  adjustcolor("lightblue", alpha = 1), 
       border = NA)
}
lines(x, y, type = "l")
lines(x, mu, col = "red", lwd = 2.5)
abline(v = x[est_cp], lwd = 2, lty = 2, col = "blue")
length(est_cp)

#### H-SMUCE 0.05 ####
fit <- stepFit(y, alpha = 0.05, confband = TRUE,
               jumpint = TRUE, family = "hsmuce") 
mu <- rep(fit$value, diff(c(fit$leftEnd,T+1)))
est_cp <- fit$leftEnd[-1]
plot(x, y, type = "l", main = TeX("H-SMUCE ($\\alpha$ = 0.05)", bold = TRUE), 
     xlab = "", ylab = "", col = "white",
     cex.main=2.5, cex.lab=2, cex.axis=1.75)
for(i in 1:length(est_cp)) { 
  rect(xleft = x[fit$leftEndLeftBound[i+1]-1]-0.5*unit, xright = x[fit$leftEndRightBound[i+1]] + 0.5*unit, 
       ybottom = par("usr")[3], ytop = par("usr")[4], 
       col =  adjustcolor("lightblue", alpha = 1), 
       border = NA)
}
lines(x, y, type = "l")
lines(x, mu, col = "red", lwd = 2.5)
abline(v = x[est_cp], lwd = 2, lty = 2, col = "blue")

#### PELT ####
fit <- cpt.meanvar(y, method = "PELT")
L_est <- length(fit@cpts) - 1
est_cp <- fit@cpts[-(L_est+1)] + 1
mu <- rep(fit@param.est$mean, diff(c(1, est_cp, T+1))) 
plot(x, y, type = "l", main = "PELT", xlab = "Time (s)", 
     ylab = "", cex.main=2.5, cex.lab=2, cex.axis=2)
lines(x, mu, col = "red", lwd = 2.5)
abline(v = x[est_cp], lwd = 2, lty = 2, col = "blue")
mtext("Conductance (nS)", side = 2, outer = TRUE, line = -1, cex = 1.5)
length(est_cp)
dev.off()

png("~/Desktop/ion_plot_2.png", width = 1250, height = 1250 / 2)
par(mfrow = c(2,1), oma = c(0,1.5,0,0), mar = c(4,4,3,2))

#### MICH Mean ####
fit_mean <- mich(y, L_auto = TRUE, max_iter = Inf, verbose = TRUE, tol = 1e-6, restart = FALSE)
cred_sets <- mich_sets(fit_mean$mean_model$pi_bar, level = 0.95)
mu <- fit_mean$mu
est_cp <- cred_sets$cp

plot(x, y, type = "l", main = paste0("MICH (L = ", fit_mean$L,")"),
     xlab = "", ylab = "", col = "white", 
     cex.main=2.5, cex.lab=2, cex.axis=2)
for(i in unlist(cred_sets$sets)) { 
  rect(xleft = x[i]-2*unit, xright = x[i] + 2*unit, 
       ybottom = par("usr")[3], ytop = par("usr")[4], 
       col =  adjustcolor("lightblue", alpha = 1), 
       border = NA)
}
# for (i in unlist(cred_sets$sets)) {
#   lines(x = c(x[i],x[i]), y = c(-1,0.075), lwd = 4, col = "blue")
# }
lines(x, y, type = "l")
lines(x, mu, col = "red", lwd = 2.5)
abline(v = x[est_cp], lwd = 2, lty = 2, col = "blue")
length(est_cp)

#### H-SMUCE 0.5 ####
fit <- stepFit(y, alpha = 0.5, confband = TRUE,
               jumpint = TRUE, family = "hsmuce") 
mu <- rep(fit$value, diff(c(fit$leftEnd,T+1)))
est_cp <- fit$leftEnd[-1]
plot(x, y, type = "l", main = TeX("H-SMUCE ($\\alpha$ = 0.5)", bold = TRUE), 
     xlab = "Time (s)", ylab = "",  
     cex.main=2.5, cex.lab=2, cex.axis=2)
for(i in 1:length(est_cp)) { 
  rect(xleft = x[fit$leftEndLeftBound[i+1]-1]-0.5*unit, xright = x[fit$leftEndRightBound[i+1]] + 0.5*unit, 
       ybottom = par("usr")[3], ytop = par("usr")[4], 
       col =  adjustcolor("lightblue", alpha = 1), 
       border = NA)
}
lines(x, y, type = "l")
lines(x, mu, col = "red", lwd = 2.5)
abline(v = x[est_cp], lwd = 2, lty = 2, col = "blue")
length(est_cp)
mtext("Conductance (nS)", side = 2, outer = TRUE, line = -1, cex = 2)

dev.off()
