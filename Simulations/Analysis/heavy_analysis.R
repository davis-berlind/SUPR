results <- readRDS("~/Desktop/SUPR/heavy_results.rds")

methods <- c("NSP", "NOT", "MOSUM BUM", "MOSUM LP", "PELT", 
             "HSMUCE (0.1)", "HSMUCE (0.5)", 
             "MICH Ora Rev", "MICH AutoFast Rev")

uq_methods <- c("NSP", "MOSUM BUM", "MOSUM LP", 
                "HSMUCE (0.1)", "HSMUCE (0.5)", 
                "MICH Ora Rev", "MICH AutoFast Rev")

png("~/Desktop/laplace_sim_dist_plot.png", width = 1100, height = 800)
results %>%
  filter(method %in% methods, family == "laplace") %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  mutate(Bias = abs(L - L_est), 
         Hausdorff = hausdorff_1 + hausdorff_2,
         FPSLE = fpsle, 
         FNSLE = fnsle, 
         ci_length = ifelse(avg_len == 0, NA, avg_len),
         Time = time,
         CCD = n_covered / n_detected) %>%  
  select(method, T, L, Bias, Hausdorff, FPSLE, FNSLE, ci_length, CCD, Time) %>% 
  pivot_longer(cols = c( Bias, Hausdorff, FPSLE, FNSLE, ci_length, CCD, Time),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(measurement = ifelse(measurement == "ci_length", "Set Length", measurement),
         alpha =  ifelse(measurement == "CCD", 0.9, NA)) %>% 
  ggplot(aes(x = as.factor(L), y = value, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  geom_hline(aes(yintercept = alpha), linetype = "dashed") +
  facet_grid(factor(measurement, levels = c("Bias", "Hausdorff", "FPSLE", "FNSLE", "Set Length", "CCD", "Time")) ~ T, 
             scales = "free", labeller = labeller(T = T.labs)) + 
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  labs(y=NULL, x = "Number of Change-Points")
dev.off()

png("~/Desktop/t_sim_dist_plot.png", width = 1100, height = 800)
results %>%
  filter(method %in% methods, family == "t") %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  mutate(Bias = abs(L - L_est), 
         Hausdorff = hausdorff_1 + hausdorff_2,
         FPSLE = fpsle, 
         FNSLE = fnsle, 
         ci_length = ifelse(avg_len == 0, NA, avg_len),
         Time = time,
         CCD = n_covered / n_detected) %>%  
  select(method, T, L, Bias, Hausdorff, FPSLE, FNSLE, ci_length, CCD, Time) %>% 
  pivot_longer(cols = c( Bias, Hausdorff, FPSLE, FNSLE, ci_length, CCD, Time),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(measurement = ifelse(measurement == "ci_length", "Set Length", measurement),
         alpha =  ifelse(measurement == "CCD", 0.9, NA)) %>% 
  ggplot(aes(x = as.factor(L), y = value, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  geom_hline(aes(yintercept = alpha), linetype = "dashed") +
  facet_grid(factor(measurement, levels = c("Bias", "Hausdorff", "FPSLE", "FNSLE", "Set Length", "CCD", "Time")) ~ T, 
             scales = "free", labeller = labeller(T = T.labs)) + 
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  labs(y=NULL, x = "Number of Change-Points")
dev.off()

png("~/Desktop/ma3_sim_dist_plot.png", width = 1100, height = 800)
results %>%
  filter(method %in% methods, family == "MA 0.3") %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  mutate(Bias = abs(L - L_est), 
         Hausdorff = hausdorff_1 + hausdorff_2,
         FPSLE = fpsle, 
         FNSLE = fnsle, 
         ci_length = ifelse(avg_len == 0, NA, avg_len),
         Time = time,
         CCD = n_covered / n_detected) %>%  
  select(method, T, L, Bias, Hausdorff, FPSLE, FNSLE, ci_length, CCD, Time) %>% 
  pivot_longer(cols = c( Bias, Hausdorff, FPSLE, FNSLE, ci_length, CCD, Time),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(measurement = ifelse(measurement == "ci_length", "Set Length", measurement),
         alpha =  ifelse(measurement == "CCD", 0.9, NA)) %>% 
  ggplot(aes(x = as.factor(L), y = value, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  geom_hline(aes(yintercept = alpha), linetype = "dashed") +
  facet_grid(factor(measurement, levels = c("Bias", "Hausdorff", "FPSLE", "FNSLE", "Set Length", "CCD", "Time")) ~ T, 
             scales = "free", labeller = labeller(T = T.labs)) + 
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  labs(y=NULL, x = "Number of Change-Points")
dev.off()

png("~/Desktop/ma7_sim_dist_plot.png", width = 1100, height = 800)
results %>%
  filter(method %in% methods, family == "MA 0.7") %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  mutate(Bias = abs(L - L_est), 
         Hausdorff = hausdorff_1 + hausdorff_2,
         FPSLE = fpsle, 
         FNSLE = fnsle, 
         ci_length = ifelse(avg_len == 0, NA, avg_len),
         Time = time,
         CCD = n_covered / n_detected) %>%  
  select(method, T, L, Bias, Hausdorff, FPSLE, FNSLE, ci_length, CCD, Time) %>% 
  pivot_longer(cols = c( Bias, Hausdorff, FPSLE, FNSLE, ci_length, CCD, Time),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(measurement = ifelse(measurement == "ci_length", "Set Length", measurement),
         alpha =  ifelse(measurement == "CCD", 0.9, NA)) %>% 
  ggplot(aes(x = as.factor(L), y = value, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  geom_hline(aes(yintercept = alpha), linetype = "dashed") +
  facet_grid(factor(measurement, levels = c("Bias", "Hausdorff", "FPSLE", "FNSLE", "Set Length", "CCD", "Time")) ~ T, 
             scales = "free", labeller = labeller(T = T.labs)) + 
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  labs(y=NULL, x = "Number of Change-Points")
dev.off()

# MA 0.3 ####
bias_star = results %>% 
  filter(method %in% methods & family == "MA 0.3") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_star = min(bias),
            hausdorff_star = min(hausdorff, na.rm = TRUE),
            fpsle_star = min(fpsle, na.rm = TRUE),
            fnsle_star = min(fnsle, na.rm = TRUE))

bias_dag = results %>% 
  filter(method %in% methods & !grepl("Ora", method) & family == "MA 0.3") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_dag = min(bias),
            hausdorff_dag = min(hausdorff, na.rm = TRUE),
            fpsle_dag = min(fpsle, na.rm = TRUE),
            fnsle_dag = min(fnsle, na.rm = TRUE))

bias_ddag = results %>% 
  filter(method %in% uq_methods & !grepl("Ora", method) & family == "MA 0.3") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_ddag = min(bias),
            hausdorff_ddag = min(hausdorff, na.rm = TRUE),
            fpsle_ddag = min(fpsle, na.rm = TRUE),
            fnsle_ddag = min(fnsle, na.rm = TRUE))

results %>% 
  filter(method %in% methods & family == "MA 0.3") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est, na.rm = TRUE),
            coverage = sum(n_covered) / sum(n_detected),
            time = mean(time, na.rm = TRUE)) %>% 
  right_join(bias_star) %>% 
  right_join(bias_dag) %>% 
  right_join(bias_ddag) %>% 
  mutate(bias_note = paste0("\\textsuperscript{", 
                            ifelse(bias == bias_star, "*", ""),
                            ifelse(bias == bias_dag, "\\dag", ""),
                            ifelse(bias == bias_ddag, "\\ddag", ""),
                            "}"),
         haus_note = paste0("\\textsuperscript{", 
                            ifelse(hausdorff == hausdorff_star, "*", ""),
                            ifelse(hausdorff == hausdorff_dag, "\\dag", ""),
                            ifelse(hausdorff == hausdorff_ddag, "\\ddag", ""),
                            "}"),
         fpsle_note = paste0("\\textsuperscript{", 
                            ifelse(fpsle == fpsle_star, "*", ""),
                            ifelse(fpsle == fpsle_dag, "\\dag", ""),
                            ifelse(fpsle == fpsle_ddag, "\\ddag", ""),
                            "}"),
         fnsle_note = paste0("\\textsuperscript{", 
                            ifelse(fnsle == fnsle_star, "*", ""),
                            ifelse(fnsle == fnsle_dag, "\\dag", ""),
                            ifelse(fnsle == fnsle_ddag, "\\ddag", ""),
                            "}"),
         ) %>% 
  mutate(ci_length = ifelse(ci_length == 0, NA, ci_length)) %>% 
  mutate_if(is.numeric, ~format(round(.,3), nsmall = 3)) %>% 
  mutate_if(is.character, ~gsub(("NaN|NA"), " ",.)) %>% 
  mutate(method = case_when(
    method == "MICH AutoFast Rev" ~ "b Auto-MICH",
    method == "MICH Ora Rev" ~ "c Ora-MICH",
    method == "HSMUCE (0.1)" ~ "d H-SMUCE (0.1)",
    method == "HSMUCE (0.5)" ~ "e H-SMUCE (0.5) ",
    method == "MOSUM BUM" ~ "f MOSUM BUM",
    method == "MOSUM LP" ~ "g MOSUM LP",
    method == "PELT" ~ "h PELT",
    method == "NOT" ~ "i NOT",
    method == "NSP" ~ "j NSP"
    )
  ) %>% 
  arrange(T, L, method) %>% 
  mutate(method = gsub("[a-z] ", "", method)) %>%
  mutate(bias = gsub(" ", "", paste0(bias_note, bias)),
         hausdorff = gsub(" ", "", paste0(haus_note, hausdorff)),
         fpsle = gsub(" ", "", paste0(fpsle_note, fpsle)), 
         fnsle = gsub(" ", "", paste0(fnsle_note, fnsle))) %>% 
  mutate_all(trimws) %>% 
  select(T, L, method, bias, hausdorff, fpsle, fnsle, ci_length, coverage, time) %>%
  write.csv("~/Desktop/SUPR/simulations/ma3_sim.csv", row.names = FALSE)

# MA 0.7 ####
bias_star = results %>% 
  filter(method %in% methods & family == "MA 0.7") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_star = min(bias),
            hausdorff_star = min(hausdorff, na.rm = TRUE),
            fpsle_star = min(fpsle, na.rm = TRUE),
            fnsle_star = min(fnsle, na.rm = TRUE))

bias_dag = results %>% 
  filter(method %in% methods & !grepl("Ora", method) & family == "MA 0.7") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_dag = min(bias),
            hausdorff_dag = min(hausdorff, na.rm = TRUE),
            fpsle_dag = min(fpsle, na.rm = TRUE),
            fnsle_dag = min(fnsle, na.rm = TRUE))

bias_ddag = results %>% 
  filter(method %in% uq_methods & !grepl("Ora", method) & family == "MA 0.7") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_ddag = min(bias),
            hausdorff_ddag = min(hausdorff, na.rm = TRUE),
            fpsle_ddag = min(fpsle, na.rm = TRUE),
            fnsle_ddag = min(fnsle, na.rm = TRUE))

results %>% 
  filter(method %in% methods & family == "MA 0.7") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est, na.rm = TRUE),
            coverage = sum(n_covered) / sum(n_detected),
            time = mean(time, na.rm = TRUE)) %>% 
  right_join(bias_star) %>% 
  right_join(bias_dag) %>% 
  right_join(bias_ddag) %>% 
  mutate(bias_note = paste0("\\textsuperscript{", 
                            ifelse(bias == bias_star, "*", ""),
                            ifelse(bias == bias_dag, "\\dag", ""),
                            ifelse(bias == bias_ddag, "\\ddag", ""),
                            "}"),
         haus_note = paste0("\\textsuperscript{", 
                            ifelse(hausdorff == hausdorff_star, "*", ""),
                            ifelse(hausdorff == hausdorff_dag, "\\dag", ""),
                            ifelse(hausdorff == hausdorff_ddag, "\\ddag", ""),
                            "}"),
         fpsle_note = paste0("\\textsuperscript{", 
                             ifelse(fpsle == fpsle_star, "*", ""),
                             ifelse(fpsle == fpsle_dag, "\\dag", ""),
                             ifelse(fpsle == fpsle_ddag, "\\ddag", ""),
                             "}"),
         fnsle_note = paste0("\\textsuperscript{", 
                             ifelse(fnsle == fnsle_star, "*", ""),
                             ifelse(fnsle == fnsle_dag, "\\dag", ""),
                             ifelse(fnsle == fnsle_ddag, "\\ddag", ""),
                             "}"),
  ) %>% 
  mutate(ci_length = ifelse(ci_length == 0, NA, ci_length)) %>% 
  mutate_if(is.numeric, ~format(round(.,3), nsmall = 3)) %>% 
  mutate_if(is.character, ~gsub(("NaN|NA"), " ",.)) %>% 
  mutate(method = case_when(
    method == "MICH AutoFast Rev" ~ "b Auto-MICH",
    method == "MICH Ora Rev" ~ "c Ora-MICH",
    method == "HSMUCE (0.1)" ~ "d H-SMUCE (0.1)",
    method == "HSMUCE (0.5)" ~ "e H-SMUCE (0.5) ",
    method == "MOSUM BUM" ~ "f MOSUM BUM",
    method == "MOSUM LP" ~ "g MOSUM LP",
    method == "PELT" ~ "h PELT",
    method == "NOT" ~ "i NOT",
    method == "NSP" ~ "j NSP"
  )
  ) %>% 
  arrange(T, L, method) %>% 
  mutate(method = gsub("[a-z] ", "", method)) %>%
  mutate(bias = gsub(" ", "", paste0(bias_note, bias)),
         hausdorff = gsub(" ", "", paste0(haus_note, hausdorff)),
         fpsle = gsub(" ", "", paste0(fpsle_note, fpsle)), 
         fnsle = gsub(" ", "", paste0(fnsle_note, fnsle))) %>% 
  mutate_all(trimws) %>% 
  select(T, L, method, bias, hausdorff, fpsle, fnsle, ci_length, coverage, time) %>%
  write.csv("~/Desktop/SUPR/simulations/ma7_sim.csv", row.names = FALSE)

# laplace ####
bias_star = results %>% 
  filter(method %in% methods & family == "laplace") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_star = min(bias),
            hausdorff_star = min(hausdorff, na.rm = TRUE),
            fpsle_star = min(fpsle, na.rm = TRUE),
            fnsle_star = min(fnsle, na.rm = TRUE))

bias_dag = results %>% 
  filter(method %in% methods & !grepl("Ora", method) & family == "laplace") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_dag = min(bias),
            hausdorff_dag = min(hausdorff, na.rm = TRUE),
            fpsle_dag = min(fpsle, na.rm = TRUE),
            fnsle_dag = min(fnsle, na.rm = TRUE))

bias_ddag = results %>% 
  filter(method %in% uq_methods & !grepl("Ora", method) & family == "laplace") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_ddag = min(bias),
            hausdorff_ddag = min(hausdorff, na.rm = TRUE),
            fpsle_ddag = min(fpsle, na.rm = TRUE),
            fnsle_ddag = min(fnsle, na.rm = TRUE))

results %>% 
  filter(method %in% methods & family == "laplace") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est, na.rm = TRUE),
            coverage = sum(n_covered) / sum(n_detected),
            time = mean(time, na.rm = TRUE)) %>% 
  right_join(bias_star) %>% 
  right_join(bias_dag) %>% 
  right_join(bias_ddag) %>% 
  mutate(bias_note = paste0("\\textsuperscript{", 
                            ifelse(bias == bias_star, "*", ""),
                            ifelse(bias == bias_dag, "\\dag", ""),
                            ifelse(bias == bias_ddag, "\\ddag", ""),
                            "}"),
         haus_note = paste0("\\textsuperscript{", 
                            ifelse(hausdorff == hausdorff_star, "*", ""),
                            ifelse(hausdorff == hausdorff_dag, "\\dag", ""),
                            ifelse(hausdorff == hausdorff_ddag, "\\ddag", ""),
                            "}"),
         fpsle_note = paste0("\\textsuperscript{", 
                             ifelse(fpsle == fpsle_star, "*", ""),
                             ifelse(fpsle == fpsle_dag, "\\dag", ""),
                             ifelse(fpsle == fpsle_ddag, "\\ddag", ""),
                             "}"),
         fnsle_note = paste0("\\textsuperscript{", 
                             ifelse(fnsle == fnsle_star, "*", ""),
                             ifelse(fnsle == fnsle_dag, "\\dag", ""),
                             ifelse(fnsle == fnsle_ddag, "\\ddag", ""),
                             "}"),
  ) %>% 
  mutate(ci_length = ifelse(ci_length == 0, NA, ci_length)) %>% 
  mutate_if(is.numeric, ~format(round(.,3), nsmall = 3)) %>% 
  mutate_if(is.character, ~gsub(("NaN|NA"), " ",.)) %>% 
  mutate(method = case_when(
    method == "MICH AutoFast Rev" ~ "b Auto-MICH",
    method == "MICH Ora Rev" ~ "c Ora-MICH",
    method == "HSMUCE (0.1)" ~ "d H-SMUCE (0.1)",
    method == "HSMUCE (0.5)" ~ "e H-SMUCE (0.5) ",
    method == "MOSUM BUM" ~ "f MOSUM BUM",
    method == "MOSUM LP" ~ "g MOSUM LP",
    method == "PELT" ~ "h PELT",
    method == "NOT" ~ "i NOT",
    method == "NSP" ~ "j NSP"
  )
  ) %>% 
  arrange(T, L, method) %>% 
  mutate(method = gsub("[a-z] ", "", method)) %>%
  mutate(bias = gsub(" ", "", paste0(bias_note, bias)),
         hausdorff = gsub(" ", "", paste0(haus_note, hausdorff)),
         fpsle = gsub(" ", "", paste0(fpsle_note, fpsle)), 
         fnsle = gsub(" ", "", paste0(fnsle_note, fnsle))) %>% 
  mutate_all(trimws) %>% 
  select(T, L, method, bias, hausdorff, fpsle, fnsle, ci_length, coverage, time) %>%
  write.csv("~/Desktop/SUPR/simulations/laplace_sim.csv", row.names = FALSE)

# t ####
bias_star = results %>% 
  filter(method %in% methods & family == "t") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_star = min(bias),
            hausdorff_star = min(hausdorff, na.rm = TRUE),
            fpsle_star = min(fpsle, na.rm = TRUE),
            fnsle_star = min(fnsle, na.rm = TRUE))

bias_dag = results %>% 
  filter(method %in% methods & !grepl("Ora", method) & family == "t") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_dag = min(bias),
            hausdorff_dag = min(hausdorff, na.rm = TRUE),
            fpsle_dag = min(fpsle, na.rm = TRUE),
            fnsle_dag = min(fnsle, na.rm = TRUE))

bias_ddag = results %>% 
  filter(method %in% uq_methods & !grepl("Ora", method) & family == "t") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L) %>% 
  summarize(bias_ddag = min(bias),
            hausdorff_ddag = min(hausdorff, na.rm = TRUE),
            fpsle_ddag = min(fpsle, na.rm = TRUE),
            fnsle_ddag = min(fnsle, na.rm = TRUE))

results %>% 
  filter(method %in% methods & family == "t") %>% 
  group_by(T, L, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est, na.rm = TRUE),
            coverage = sum(n_covered) / sum(n_detected),
            time = mean(time, na.rm = TRUE)) %>% 
  right_join(bias_star) %>% 
  right_join(bias_dag) %>% 
  right_join(bias_ddag) %>% 
  mutate(bias_note = paste0("\\textsuperscript{", 
                            ifelse(bias == bias_star, "*", ""),
                            ifelse(bias == bias_dag, "\\dag", ""),
                            ifelse(bias == bias_ddag, "\\ddag", ""),
                            "}"),
         haus_note = paste0("\\textsuperscript{", 
                            ifelse(hausdorff == hausdorff_star, "*", ""),
                            ifelse(hausdorff == hausdorff_dag, "\\dag", ""),
                            ifelse(hausdorff == hausdorff_ddag, "\\ddag", ""),
                            "}"),
         fpsle_note = paste0("\\textsuperscript{", 
                             ifelse(fpsle == fpsle_star, "*", ""),
                             ifelse(fpsle == fpsle_dag, "\\dag", ""),
                             ifelse(fpsle == fpsle_ddag, "\\ddag", ""),
                             "}"),
         fnsle_note = paste0("\\textsuperscript{", 
                             ifelse(fnsle == fnsle_star, "*", ""),
                             ifelse(fnsle == fnsle_dag, "\\dag", ""),
                             ifelse(fnsle == fnsle_ddag, "\\ddag", ""),
                             "}"),
  ) %>% 
  mutate(ci_length = ifelse(ci_length == 0, NA, ci_length)) %>% 
  mutate_if(is.numeric, ~format(round(.,3), nsmall = 3)) %>% 
  mutate_if(is.character, ~gsub(("NaN|NA"), " ",.)) %>% 
  mutate(method = case_when(
    method == "MICH AutoFast Rev" ~ "b Auto-MICH",
    method == "MICH Ora Rev" ~ "c Ora-MICH",
    method == "HSMUCE (0.1)" ~ "d H-SMUCE (0.1)",
    method == "HSMUCE (0.5)" ~ "e H-SMUCE (0.5) ",
    method == "MOSUM BUM" ~ "f MOSUM BUM",
    method == "MOSUM LP" ~ "g MOSUM LP",
    method == "PELT" ~ "h PELT",
    method == "NOT" ~ "i NOT",
    method == "NSP" ~ "j NSP"
  )
  ) %>% 
  arrange(T, L, method) %>% 
  mutate(method = gsub("[a-z] ", "", method)) %>%
  mutate(bias = gsub(" ", "", paste0(bias_note, bias)),
         hausdorff = gsub(" ", "", paste0(haus_note, hausdorff)),
         fpsle = gsub(" ", "", paste0(fpsle_note, fpsle)), 
         fnsle = gsub(" ", "", paste0(fnsle_note, fnsle))) %>% 
  mutate_all(trimws) %>% 
  select(T, L, method, bias, hausdorff, fpsle, fnsle, ci_length, coverage, time) %>%
  write.csv("~/Desktop/SUPR/simulations/t_sim.csv", row.names = FALSE)
