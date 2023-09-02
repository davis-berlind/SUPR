#### plot functions #### 

plot_mich <- function(mich_fit, level = 0.95) {
  alph <- (1 - level) / 2
  q <- -qnorm(alph)
  plot(mich_fit$y, ylim = c(min(mich_fit$y - q / sqrt(mich_fit$lambda)), max(mich_fit$y + q / sqrt(mich_fit$lambda))),
       ylab = "")
  lines(mich_fit$mu, col = "red", lwd = 3)
}
