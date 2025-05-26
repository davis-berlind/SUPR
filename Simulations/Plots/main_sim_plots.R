library(tidyverse)
library(ggplot2)
library(paletteer)

results <- readRDS("~/Desktop/SUPR/results.rds")

methods <- c("NSP", "MOSUM LP", "PELT", 
             "HSMUCE (0.1)", "HSMUCE (0.5)", 
             "MICH Ora Rev", "MICH AutoFast Rev")

uq_methods <- c("NSP", "MOSUM BUM", "MOSUM LP", 
                "HSMUCE (0.1)", "HSMUCE (0.5)", 
                "MICH Ora Rev", "MICH Auto Rev", "MICH AutoFast Rev")

T.labs <- paste0("T = ", c(100, 500, 1000)) 
names(T.labs) <- c(100, 500, 1000)

# slide plots ####
png("~/Desktop/sim_plot_slides.png", width = 1425, height = 1425/2)
results %>% 
  filter(method %in% c("PELT",  "HSMUCE (0.1)", "HSMUCE (0.5)", 
                       "MICH Ora Rev", "MICH AutoFast Rev")) %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-10 | is.na(delta))) %>% 
  group_by(T, L, method) %>% 
  summarize(Bias = mean(abs(L - L_est)),
            FPSLE = mean(fpsle),
            FNSLE = mean(fnsle),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est)) %>%
  mutate(ci_length = ifelse(ci_length == 0, NA, ci_length)) %>% 
  select(method, T, L, Bias, FPSLE, FNSLE, ci_length) %>% 
  pivot_longer(cols = c(Bias, FPSLE, FNSLE, ci_length),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(measurement = ifelse(measurement == "ci_length", "CI Length", measurement)) %>% 
  ggplot(aes(x = L, y = value, color = method, group = method)) + 
  facet_grid(factor(measurement, levels = c("Bias", "FPSLE", "FNSLE", "CI Length")) ~ T, 
             scales = "free", labeller = labeller(T = T.labs)) + 
  geom_line(size = 1.5, alpha = 0.8, linetype = "solid") +
  geom_point(aes(shape = method, fill = method), alpha = 0.5, size = 3) + 
  scale_shape_manual(values=c(25:21)) +
  # scale_colour_paletteer_d("ggthemes::Red_Blue_Brown") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
labs(y=NULL)
dev.off()

png("~/Desktop/sim_plot_med_slides.png", width = 1425, height = 1425/2)
results %>% 
  filter(method %in% c("PELT",  "HSMUCE (0.1)", "HSMUCE (0.5)", 
                       "MICH Ora Rev", "MICH AutoFast Rev")) %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-10 | is.na(delta))) %>% 
  group_by(T, L, method) %>% 
  summarize(Bias = median(abs(L - L_est)),
            FPSLE = median(fpsle),
            FNSLE = median(fnsle),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est)) %>%
  mutate(ci_length = ifelse(ci_length == 0, NA, ci_length)) %>% 
  select(method, T, L, Bias, FPSLE, FNSLE, ci_length) %>% 
  pivot_longer(cols = c(Bias, FPSLE, FNSLE, ci_length),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(measurement = ifelse(measurement == "ci_length", "CI Length", measurement)) %>% 
  ggplot(aes(x = L, y = value, color = method, group = method)) + 
  facet_grid(factor(measurement, levels = c("Bias", "FPSLE", "FNSLE", "CI Length")) ~ T, 
             scales = "free", labeller = labeller(T = T.labs)) + 
  geom_line(size = 1.5, alpha = 0.8, linetype = "solid") +
  geom_point(aes(shape = method, fill = method), alpha = 0.5, size = 3) + 
  scale_shape_manual(values=c(25:21)) +
  # scale_colour_paletteer_d("ggthemes::Red_Blue_Brown") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  labs(y=NULL)
dev.off()

png("~/Desktop/sim_plot.png", width = 1100, height = 600)
results %>% 
  filter(method %in% methods) %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
  group_by(T, L, method) %>% 
  summarize(Bias = mean(abs(L - L_est)),
            # Hausdorff = mean(hausdorff_1 + hausdorff_2),
            FPSLE = mean(fpsle),
            FNSLE = mean(fnsle),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est)) %>%
  mutate(ci_length = ifelse(ci_length == 0, NA, ci_length)) %>% 
  select(method, T, L, Bias, FPSLE, FNSLE, ci_length) %>% 
  pivot_longer(cols = c(Bias, FPSLE, FNSLE, ci_length),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(measurement = ifelse(measurement == "ci_length", "Set Length", measurement)) %>% 
  ggplot(aes(x = L, y = value, color = method, group = method)) + 
  facet_grid(factor(measurement, levels = c("Bias", "FPSLE", "FNSLE", "Set Length")) ~ T, 
             scales = "free", labeller = labeller(T = T.labs)) + 
  geom_line(size = 1.5, alpha = 1, linetype = "solid") +
  geom_point(aes(shape = method, fill = method), alpha = 0.5, size = 2.5) + 
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("MetBrewer::Austria") +
  scale_fill_paletteer_d("MetBrewer::Austria") +
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

png("~/Desktop/sim_med_plot.png", width = 1100, height = 600)
results %>% 
  filter(method %in% methods) %>% 
  mutate(method = gsub("MICH AutoFast Rev", "Auto-MICH", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "Ora-MICH", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
  group_by(T, L, method) %>% 
  summarize(Bias = median(abs(L - L_est)),
            # Hausdorff = mean(hausdorff_1 + hausdorff_2),
            FPSLE = median(fpsle),
            FNSLE = median(fnsle),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est)) %>%
  mutate(ci_length = ifelse(ci_length == 0, NA, ci_length)) %>% 
  select(method, T, L, Bias, FPSLE, FNSLE, ci_length) %>% 
  pivot_longer(cols = c(Bias, FPSLE, FNSLE, ci_length),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(measurement = ifelse(measurement == "ci_length", "Set Length", measurement)) %>% 
  ggplot(aes(x = as.factor(L), y = value, color = method, fill = method)) + 
  geom_boxplot() + 
  facet_grid(factor(measurement, levels = c("Bias", "FPSLE", "FNSLE", "Set Length")) ~ T, 
             scales = "free", labeller = labeller(T = T.labs)) + 
  geom_line(size = 1.5, alpha = 1, linetype = "solid") +
  geom_point(aes(shape = method, fill = method), alpha = 0.5, size = 2.5) +
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("MetBrewer::Austria") +
  scale_fill_paletteer_d("MetBrewer::Austria") +
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

methods <- c("PELT", 
             "HSMUCE (0.1)", "HSMUCE (0.5)", 
             "MICH Ora Rev", "MICH AutoFast Rev")

png("~/Desktop/sim_dist_plot.png", width = 1100, height = 650)
results %>% 
  filter(method %in% methods) %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
  mutate(Bias = abs(L - L_est), 
         FPSLE = fpsle, 
         FNSLE = fnsle, 
         `Set Length` = ifelse(avg_len == 0, NA, avg_len),
         `Time (s)` = time,
         CCD = n_covered / n_detected) %>%  
  select(method, T, L, Bias, FPSLE, FNSLE, `Set Length`, CCD, `Time (s)`) %>% 
  pivot_longer(cols = c( Bias, FPSLE, FNSLE, `Set Length`, CCD, `Time (s)`),
               names_to = "measurement", 
               values_to = "value") %>% 
  mutate(alpha =  ifelse(measurement == "CCD", 0.9, NA)) %>% 
  ggplot(aes(x = as.factor(L), y = value, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5, stroke = 1.5) + 
  geom_hline(aes(yintercept = alpha), linetype = "dashed") +
  facet_grid(factor(measurement, levels = c("Bias", "FPSLE", "FNSLE", "Set Length", "CCD", "Time (s)")) ~ T, 
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

methods <- c("NSP", "MOSUM BUM", "MOSUM LP", "NOT", "PELT", 
             "HSMUCE (0.1)", "HSMUCE (0.5)", 
             "MICH Ora Rev", "MICH AutoFast Rev")

png("~/Desktop/low_full_sim_dist_plot.png", width = 1100, height = 800)
results %>% 
  filter(method %in% methods) %>% 
  filter((T == 100 & min_space == 15) | (T == 500 & min_space == 15) | (T == 1000 & min_space == 30)) %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
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

png("~/Desktop/high_full_sim_dist_plot.png", width = 1100, height = 800)
results %>% 
  filter(method %in% methods) %>% 
  filter((T == 500 & min_space == 30) | (T == 1000 & min_space == 50)) %>% 
  mutate(method = gsub("MICH AutoFast Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("HSMUCE", "H-SMUCE", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
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
         alpha = ifelse(measurement == "CCD", 0.9, NA)) %>% 
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

