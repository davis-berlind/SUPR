library(latex2exp)

x = seq(0,10,length.out=10000)
plot(x, x -log(x) -1, type = "l", lwd = 3,
     ylab = "Signal Strength",
     xlab = TeX("$s_0^2$"),
     main = "Variance Signal Control Functions")
x = seq(0,10,length.out=500)
lines(x, 1 / x + log(x) - 1, col = "blue", lty = 2, lwd = 3)
legend("topright", legend = c( TeX("$f_1$"),  TeX("$f_2$")), lwd = 2, 
       col = c("black", "blue"),lty = c(1,2) , cex = 2)
