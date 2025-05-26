set.seed(11)
T = 50
B_r = 0
d = 1
N = 5000
y = sapply(1:N, function(t) rnorm(T))

log_pi_mean_prior = log(mean_prior(T, B_r, d))
pi_mean = sapply(1:N,
                 function(t) mean_scp(y[,t], lambda = rep(1,T), 
                                      tau = 0, log_pi = rep(-log(T),T), 
                                      B_r = 0)$pi_bar)
pi_mean_wtd = sapply(1:N, 
                     function(t) mean_scp(y[,t], lambda = rep(1,T),
                                          tau = 0, log_pi = log_pi_mean_prior, 
                                          B_r = 0)$pi_bar)

log_pi_var_prior = log(var_prior(T, B_r))
pi_var = sapply(1:N, 
                function(t) var_scp(y[,t], tau = rep(1,T), u = 0, v = rep(0, T),
                                    log_pi = rep(-log(T),T), B_r = 0)$pi_bar)
pi_var_wtd = sapply(1:N, 
                    function(t) var_scp(y[,t], tau = rep(1,T), 
                                        u = 0, v = rep(0, T), 
                                        log_pi = log_pi_var_prior,
                                        B_r = 0)$pi_bar)

log_pi_meanvar_prior = log(meanvar_prior(T, B_r))
pi_meanvar = sapply(1:N, 
                    function(t) meanvar_scp(y[,t], lambda = rep(1,T), 
                                            tau = 0, u = 0, v = rep(0.001, T), 
                                            log_pi = rep(-log(T),T), 
                                            B_r = 0)$pi_bar)
pi_meanvar_wtd = sapply(1:N, 
                        function(t) meanvar_scp(y[,t], lambda = rep(1,T), 
                                                tau = 0, u = 0, v = rep(0.001, T), 
                                                log_pi = log_pi_meanvar_prior, 
                                                B_r = 0)$pi_bar)

par(mfrow = c(3,1), mar = c(2,2,1,0))
mean_df = data.frame(prior = pi_mean_prior,
                     uniform = apply(pi_mean, 1, mean),
                     weighted = apply(pi_mean_wtd, 1, mean))
barplot(t(mean_df), beside = TRUE, col = c("red", "grey", "black"),
        names.arg = 1:T, cex.names = 0.5,
        main = "Mean-SCP Model")
abline(h = 1 / T, lty = 3)
legend("topleft", legend = c("Weighted Prior",
                             "Posterior - Uniform Prior", 
                             "Posterior - Weighted Prior"), 
       fill = c("red", "grey","black"), cex = 1, 
       bty = "n")

var_df = data.frame(prior = pi_var_prior,
                    uniform = apply(pi_var, 1, mean),
                    weighted = apply(pi_var_wtd, 1, mean))
barplot(t(var_df), beside = TRUE, col = c("red", "grey", "black"),
        names.arg = 1:T, cex.names = 0.5,
        main = "Var-SCP Model")
abline(h = 1 / T, lty = 3)

meanvar_df = data.frame(prior = pi_meanvar_prior, 
                        uniform = apply(pi_meanvar, 1, mean),
                        weighted = apply(pi_meanvar_wtd, 1, mean))
barplot(t(meanvar_df), beside = TRUE, col = c("red", "grey", "black"),
        names.arg = 1:T, cex.names = 0.5,
        main = "MeanVar-SCP Model")
abline(h = 1 / T, lty = 3)


