results <- readRDS("~/Desktop/SUPR/simulations/results.rds")

results %>% 
  filter(method %in% "MICH Ora Rev") %>% 
  filter(delta == 0.5, tol == 1e-7) %>% 
  group_by(T, L, min_space, prior) %>% 
  summarize("|J - J_hat|" = mean(abs(L - L_est)),
            # bias_se = sd(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            # hausdorff_se = sd(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            # fpsle_se = sd(fpsle, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            # fnsle_se = sd(fnsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est, na.rm = TRUE),
            coverage = sum(n_covered) / sum(n_detected),
            time = mean(time, na.rm = TRUE)) %>% 
  mutate_if(is.numeric, ~format(round(.,3), nsmall = 3)) %>% 
  arrange(T, L, min_space, prior) %>% 
  write.csv("~/Desktop/SUPR/simulations/ora_prior_sa.csv", row.names = FALSE)

results %>% 
  filter(method %in% "MICH AutoFast Rev") %>% 
  filter(delta == 0.5, tol == 1e-7) %>% 
  group_by(T, L, min_space, prior) %>% 
  summarize("|J - J_hat|" = mean(abs(L - L_est)),
            # bias_se = sd(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            # hausdorff_se = sd(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            # fpsle_se = sd(fpsle, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            # fnsle_se = sd(fnsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est, na.rm = TRUE),
            coverage = sum(n_covered) / sum(n_detected),
            time = mean(time, na.rm = TRUE)) %>% 
  mutate_if(is.numeric, ~format(round(.,3), nsmall = 3)) %>% 
  arrange(T, L, min_space, prior) %>% 
  write.csv("~/Desktop/SUPR/simulations/auto_prior_sa.csv", row.names = FALSE)

