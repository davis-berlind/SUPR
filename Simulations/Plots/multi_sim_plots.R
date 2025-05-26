library(tidyverse)
library(ggplot2)
library(paletteer)

results <- readRDS("~/Desktop/SUPR/multi_results.rds")

methods <- c("Inspect", "L2HDC", "ecp", "MICH Ora rev", "MICH Auto rev")

results <- results %>% 
  filter(method %in% methods) %>% 
  mutate(method = gsub("MICH Auto rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("ecp", "E-Divisive", method)) %>% 
  mutate(method = gsub("L2HDC", "L2-HD", method)) %>% 
  mutate(method = gsub("MICH Ora rev", "MICH-Ora", method)) 

# main plot ####
bias_plot <- results %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  mutate(Bias = abs(L - L_est)) %>%  
  select(method, L, d, p, Bias) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = Bias, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Bias", y = NULL, x = NULL)

time_plot <- results %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  select(method, L, d, p, time) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = time, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Time (s)", y = NULL, x = NULL)

fpsle_plot <- results %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  select(method, L, d, p, fpsle) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = fpsle, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "FPSLE", y = NULL, x = NULL)

fnsle_plot <- results %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  select(method, L, d, p, fnsle) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = fnsle, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "FNSLE", y = NULL, x = NULL)

ccd_plot <- results %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  mutate(CCD = n_covered / n_detected, alpha = 0.9) %>% 
  select(method, L, d, p, CCD, alpha) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = CCD, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_fill_manual(values = c("#117733FF","#332288FF")) + 
  scale_color_manual(values = c("#117733FF","#332288FF")) +   
  geom_hline(aes(yintercept = alpha), linetype = "dashed") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "CCD", y = NULL, x = NULL)

length_plot <- results %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  mutate(set_length = ifelse(avg_len == 0, NA, avg_len)) %>% 
  select(method, L, d, p, set_length) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = set_length, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_fill_manual(values = c("#117733FF","#332288FF")) + 
  scale_color_manual(values = c("#117733FF","#332288FF")) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Set Length", y = NULL, x = NULL)

png("~/Desktop/multi_sim_dist_plot.png", width = 1500, height = 1800)
ggpubr::ggarrange(bias_plot, time_plot, 
                  fpsle_plot, fnsle_plot,
                  length_plot, ccd_plot,
                  nrow=3, ncol = 2, 
                  common.legend = TRUE, legend="bottom")

dev.off()

results %>% 
  filter(method %in% methods) %>% 
  mutate(method = gsub("MICH Auto Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("ecp", "E-Divisive", method)) %>% 
  mutate(method = gsub("L2HDC", "L2-HD", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  mutate(Bias = abs(L - L_est), 
         Hausdorff = hausdorff_1 + hausdorff_2,
         FPSLE = fpsle, 
         FNSLE = fnsle, 
         `Set Length` = ifelse(avg_len == 0, NA, avg_len),
         `Time (s)` = time,
         CCD = n_covered / n_detected) %>%  
  select(method, L, d, p, Bias) %>% 
  mutate(alpha = ifelse(measurement == "CCD", 0.9, NA)) %>% 
  ggplot(aes(x = as.factor(p), y = value, color = method, fill = method)) + 
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

results %>% 
  filter(method %in% methods) %>% 
  mutate(method = gsub("MICH Auto Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("ecp", "E-Divisive", method)) %>% 
  mutate(method = gsub("L2HDC", "L2-HD", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  mutate(Bias = abs(L - L_est), 
         Hausdorff = hausdorff_1 + hausdorff_2,
         FPSLE = fpsle, 
         FNSLE = fnsle, 
         `Set Length` = ifelse(avg_len == 0, NA, avg_len),
         `Time (s)` = time,
         CCD = n_covered / n_detected) %>%  
  select(method, L, d, p, Bias) %>% 
  mutate(alpha = ifelse(measurement == "CCD", 0.9, NA)) %>% 
  ggplot(aes(x = as.factor(p), y = value, color = method, fill = method)) + 
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

results %>% 
  filter(method %in% methods) %>% 
  mutate(method = gsub("MICH Auto Rev", "MICH-Auto", method)) %>% 
  mutate(method = gsub("ecp", "E-Divisive", method)) %>% 
  mutate(method = gsub("L2HDC", "L2-HD", method)) %>% 
  mutate(method = gsub("MICH Ora Rev", "MICH-Ora", method)) %>% 
  filter(rho == 0, adapt == FALSE) %>% 
  mutate(Bias = abs(L - L_est), 
         Hausdorff = hausdorff_1 + hausdorff_2,
         FPSLE = fpsle, 
         FNSLE = fnsle, 
         `Set Length` = ifelse(avg_len == 0, NA, avg_len),
         `Time (s)` = time,
         CCD = n_covered / n_detected) %>%  
  select(method, L, d, p, Bias) %>% 
  mutate(alpha = ifelse(measurement == "CCD", 0.9, NA)) %>% 
  ggplot(aes(x = as.factor(p), y = value, color = method, fill = method)) + 
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

# spatial correlation plot ####
bias_plot <- results %>% 
  filter(rho == 0.7, adapt == FALSE) %>% 
  mutate(Bias = abs(L - L_est)) %>%  
  select(method, L, d, p, Bias) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = Bias, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Bias", y = NULL, x = NULL)

time_plot <- results %>% 
  filter(rho == 0.7, adapt == FALSE) %>% 
  select(method, L, d, p, time) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = time, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Time (s)", y = NULL, x = NULL)

fpsle_plot <- results %>% 
  filter(rho == 0.7, adapt == FALSE) %>% 
  select(method, L, d, p, fpsle) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = fpsle, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "FPSLE", y = NULL, x = NULL)

fnsle_plot <- results %>% 
  filter(rho == 0.7, adapt == FALSE) %>% 
  select(method, L, d, p, fnsle) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = fnsle, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "FNSLE", y = NULL, x = NULL)

ccd_plot <- results %>% 
  filter(rho == 0.7, adapt == FALSE) %>% 
  mutate(CCD = n_covered / n_detected, alpha = 0.9) %>% 
  select(method, L, d, p, CCD, alpha) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = CCD, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_fill_manual(values = c("#117733FF","#332288FF")) + 
  scale_color_manual(values = c("#117733FF","#332288FF")) +   
  geom_hline(aes(yintercept = alpha), linetype = "dashed") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "CCD", y = NULL, x = NULL)

length_plot <- results %>% 
  filter(rho == 0.7, adapt == FALSE) %>% 
  mutate(set_length = ifelse(avg_len == 0, NA, avg_len)) %>% 
  select(method, L, d, p, set_length) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = set_length, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_fill_manual(values = c("#117733FF","#332288FF")) + 
  scale_color_manual(values = c("#117733FF","#332288FF")) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Set Length", y = NULL, x = NULL)

png("~/Desktop/corr_multi_sim_dist_plot.png", width = 1500, height = 1800)
ggpubr::ggarrange(bias_plot, time_plot, 
                  fpsle_plot, fnsle_plot,
                  length_plot, ccd_plot,
                  nrow=3, ncol = 2, 
                  common.legend = TRUE, legend="bottom")

dev.off()


# adapt plot ####
bias_plot <- results %>% 
  filter(rho == 0, adapt == TRUE) %>% 
  mutate(Bias = abs(L - L_est)) %>%  
  select(method, L, d, p, Bias) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = Bias, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_shape_manual(values=c(25:21,25:24)) +  
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Bias", y = NULL, x = NULL)

time_plot <- results %>% 
  filter(rho == 0, adapt == TRUE) %>% 
  select(method, L, d, p, time) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = time, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Time (s)", y = NULL, x = NULL)

fpsle_plot <- results %>% 
  filter(rho == 0, adapt == TRUE) %>% 
  select(method, L, d, p, fpsle) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = fpsle, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "FPSLE", y = NULL, x = NULL)

fnsle_plot <- results %>% 
  filter(rho == 0, adapt == TRUE) %>% 
  select(method, L, d, p, fnsle) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = fnsle, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_color_paletteer_d("rcartocolor::Safe") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "FNSLE", y = NULL, x = NULL)

ccd_plot <- results %>% 
  filter(rho == 0, adapt == TRUE) %>% 
  mutate(CCD = n_covered / n_detected, alpha = 0.9) %>% 
  select(method, L, d, p, CCD, alpha) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = CCD, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_fill_manual(values = c("#117733FF","#332288FF")) + 
  scale_color_manual(values = c("#117733FF","#332288FF")) +   
  geom_hline(aes(yintercept = alpha), linetype = "dashed") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "CCD", y = NULL, x = NULL)

length_plot <- results %>% 
  filter(rho == 0, adapt == TRUE) %>% 
  mutate(set_length = ifelse(avg_len == 0, NA, avg_len)) %>% 
  select(method, L, d, p, set_length) %>% 
  ggplot(aes(x = factor(p, labels = paste0("p = ", c(0.1,0.5,1))), y = set_length, color = method, fill = method)) + 
  geom_boxplot(outliers = FALSE, alpha = 0.2, fatten = 3) +
  stat_summary(fun.y="mean", geom="point", position = position_dodge(width = .75), shape = 5, size = 2.5) + 
  facet_grid(factor(d, labels = paste0("d = ", c(10,50,100))) ~ factor(L, labels = paste0("L = ", c(5,10,20))) , scales = "free") + 
  scale_fill_manual(values = c("#117733FF","#332288FF")) + 
  scale_color_manual(values = c("#117733FF","#332288FF")) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=24, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=24)) +
  labs(title = "Set Length", y = NULL, x = NULL)

png("~/Desktop/adapt_multi_sim_dist_plot.png", width = 1500, height = 1800)
ggpubr::ggarrange(bias_plot, time_plot, 
                  fpsle_plot, fnsle_plot,
                  length_plot, ccd_plot,
                  nrow=3, ncol = 2, 
                  common.legend = TRUE, legend="bottom")

dev.off()

