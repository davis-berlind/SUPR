library(ggplot2)
library(paletteer)
library(tidyverse)
library(InspectChangepoint)
library(L2hdchange)
library(ecp)

standardize <- function(df) {
  df <- as.matrix(df)
  varcov <- var(apply(df, 2, diff)) / 2
  r <- eigen(varcov)
  df %*% r$vectors %*% diag(1/sqrt(r$values)) %*% t(r$vectors)
}

well_log <- read.csv("~/Desktop/facies_data.csv") 

well_log %>% 
  group_by(Well.Name) %>% 
  summarise(changes = sum(diff(Facies) != 0))

wells <- unique(well_log$Well.Name)

well <- well_log %>% 
  filter(grepl("(B|C)", Formation)) %>% 
  filter(Well.Name == "SHANKLE") %>% 
  select(Facies, Depth, GR, ILD_log10, DeltaPHI, PHIND, PE, NM_M) %>% 
  rename(MnM = NM_M, Rt = ILD_log10, DeltaPhi = DeltaPHI, AvgPhi = PHIND)

true_cp <- which(diff(well$Facies) != 0) + 1
T <- nrow(well)

well <- well %>% mutate(Facies = case_when(Facies == 1 ~ "Nonmarine Sandstone",
                                           Facies == 2 ~ "Nonmarine Coarse Siltstone",
                                           Facies == 3 ~ "Nonmarine Fine Siltstone",
                                           Facies == 4 ~ "Marine Siltstone/Shale",
                                           Facies == 5 ~ "Mudstone",
                                           Facies == 6 ~ "Wackestone",
                                           Facies == 7 ~ "Dolomite",
                                           Facies == 8 ~ "Packstone-Grainstone",
                                           Facies == 9 ~ "Phylloid-Algal Bafflestone"))

well_pivot <- well[-nrow(well), ] %>%
  pivot_longer(cols = c(GR, Rt, DeltaPhi, AvgPhi, PE, MnM),
               names_to = "Measurement", 
               values_to = "Value") 
well_pivot2<- well[-1, ] %>%
  pivot_longer(cols = c(GR, Rt, DeltaPhi, AvgPhi, PE, MnM),
               names_to = "Measurement", 
               values_to = "Value") 

png("~/Desktop/well_log.png", width = 1600, height = (6.5 / 10) * 1300)
well_pivot %>%
  mutate(Value_end = well_pivot2$Value,
         Depth_end = well_pivot2$Depth) %>% 
  ggplot() + 
  geom_segment(aes(x = Depth, xend = Depth_end, y = Value, yend = Value_end, color = Facies), size = 4.5) +
  geom_point(aes(x = Depth, y = Value, color = Facies), size = 4) +
  scale_colour_paletteer_d("tvthemes::Bismuth") +
  facet_grid(rows = vars(Measurement), scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position="bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20)) +
  labs(title = "Lithology of Shankle Oil Well",
       y=NULL, x = "Depth (ft)")
dev.off()

# fit MICH ####
fit <- mich(well[,-c(1,2)], L_auto = TRUE, tol = 1e-10, restart = FALSE, verbose = TRUE)
fit_rev <- mich(well[,-c(1,2)], L_auto = TRUE, tol = 1e-10, restart = FALSE, reverse = TRUE, verbose = TRUE)

if (max(fit_rev$elbo) > max(fit$elbo)) fit <- fit_rev

est_cp <- mich_sets(fit_fast_rev$pi_bar_l, level = 0.99)$cp
sets <- mich_sets(fit_fast_rev$pi_bar_l, level = 0.99)$sets

sum(apply(abs(outer(unlist(sets), true_cp, `-`)), 2, min) <= 1)

length(est_cp[sapply(sets, function(set) min(apply(abs(outer(set, true_cp, `-`)),1,min))) <= 0]) 
length(est_cp[sapply(sets, function(set) min(apply(abs(outer(set, true_cp, `-`)),1,min))) <= 1]) 

png("~/Desktop/mich_well_log.png", width = 1300, height = (6.5 / 10) * 1300)
well_pivot %>%
  mutate(Value_end = well_pivot2$Value,
         Depth_end = well_pivot2$Depth) %>% 
  ggplot() + 
  geom_vline(xintercept = well$Depth[unlist(sets)], color = "lightblue", alpha = 0.6, size = 2.5) +
  geom_segment(aes(x = Depth, xend = Depth_end, y = Value, yend = Value_end, color = Facies), size = 4.5) +
  geom_point(aes(x = Depth, y = Value, color = Facies), size = 4) +
  scale_colour_paletteer_d("tvthemes::Bismuth") +
  geom_vline(xintercept = well$Depth[est_cp[sapply(sets, function(set) min(apply(abs(outer(set, true_cp, `-`)),1,min))) > 1]], 
             linetype = "dashed", linewidth = 1.1, color = "black", alpha = 0.6) +
  geom_vline(xintercept =  well$Depth[est_cp[sapply(sets, function(set) min(apply(abs(outer(set, true_cp, `-`)),1,min))) <= 1]], 
             linetype = "solid", linewidth = 1.25, color = "black", alpha = 0.6) +
  facet_grid(rows = vars(Measurement), scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position="bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20)) +
  labs(title = "Lithology of Shankle Oil Well",
       y=NULL, x = "Depth (ft)")

dev.off()

# fit inspect ####
inspect_fit <- inspect(t(as.matrix(well[,-c(1,2,8)])))
est_cp = inspect_fit$changepoints[,1]
sum(apply(abs(outer(est_cp, true_cp, `-`)), 2, min) <= 1)

png("~/Desktop/inspect_well_log.png", width = 1300, height = (6.5 / 10) * 1300)
well_pivot %>%
  mutate(Value_end = well_pivot2$Value,
         Depth_end = well_pivot2$Depth) %>% 
  ggplot() + 
  geom_segment(aes(x = Depth, xend = Depth_end, y = Value, yend = Value_end, color = Facies), size = 4.5) +
  geom_point(aes(x = Depth, y = Value, color = Facies), size = 4) +
  scale_colour_paletteer_d("tvthemes::Bismuth") +
  geom_vline(xintercept = well$Depth[est_cp[apply(abs(outer(est_cp, true_cp, `-`)),1,min) > 1]], 
             linetype = "dashed", linewidth = 1.1, color = "black", alpha = 0.6) +
  geom_vline(xintercept =  well$Depth[est_cp[apply(abs(outer(est_cp, true_cp, `-`)),1,min) <= 1]], 
             linetype = "solid", linewidth = 1.25, color = "black", alpha = 0.6) +
  facet_grid(rows = vars(Measurement), scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position="bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20)) +
  labs(title = "Lithology of Shankle Oil Well",
       y=NULL, x = "Depth (ft)")
dev.off()

# fit l2hdchange ####
ts_l2_fit <- ts_hdchange(t(as.matrix(well[,-c(1,2)])))
l2_fit <- hdchange(ts_l2_fit)
est_cp = l2_fit$time_stamps

png("~/Desktop/l2hdc_well_log.png", width = 1300, height = (6.5 / 10) * 1300)
well_pivot %>%
  mutate(Value_end = well_pivot2$Value,
         Depth_end = well_pivot2$Depth) %>% 
  ggplot() + 
  geom_segment(aes(x = Depth, xend = Depth_end, y = Value, yend = Value_end, color = Facies), size = 4.5) +
  geom_point(aes(x = Depth, y = Value, color = Facies), size = 4) +
  scale_colour_paletteer_d("tvthemes::Bismuth") +
  geom_vline(xintercept = well$Depth[est_cp[apply(abs(outer(est_cp, true_cp, `-`)),1,min) > 1]], 
             linetype = "dashed", linewidth = 1.1, color = "black", alpha = 0.6) +
  geom_vline(xintercept =  well$Depth[est_cp[apply(abs(outer(est_cp, true_cp, `-`)),1,min) <= 1]], 
             linetype = "solid", linewidth = 1.25, color = "black", alpha = 0.6) +
  facet_grid(rows = vars(Measurement), scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position="bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20)) +
  labs(title = "Lithology of Shankle Oil Well",
       y=NULL, x = "Depth (ft)")
dev.off()

# fit eco ####
inspect_fit <- e.divisive(as.matrix(well[,-c(1,2)]), sig.lvl = 0.01, 
                          min.size = min_space, alpha = 2,
                          R=499)
est_cp <- inspect_fit$estimates[-c(1, length(inspect_fit$estimates))]
sum(apply(abs(outer(est_cp, true_cp, `-`)), 2, min) <= 1)

png("~/Desktop/ecp_well_log.png", width = 1300, height = (6.5 / 10) * 1300)
well_pivot %>%
  mutate(Value_end = well_pivot2$Value,
         Depth_end = well_pivot2$Depth) %>% 
  ggplot() + 
  geom_segment(aes(x = Depth, xend = Depth_end, y = Value, yend = Value_end, color = Facies), size = 4.5) +
  geom_point(aes(x = Depth, y = Value, color = Facies), size = 4) +
  scale_colour_paletteer_d("tvthemes::Bismuth") +
  geom_vline(xintercept = well$Depth[est_cp[apply(abs(outer(est_cp, true_cp, `-`)),1,min) > 1]], 
             linetype = "dashed", linewidth = 1.1, color = "black", alpha = 0.6) +
  geom_vline(xintercept =  well$Depth[est_cp[apply(abs(outer(est_cp, true_cp, `-`)),1,min) <= 1]], 
             linetype = "solid", linewidth = 1.25, color = "black", alpha = 0.6) +
  facet_grid(rows = vars(Measurement), scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position="bottom", 
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20)) +
  labs(title = "Lithology of Shankle Oil Well",
       y=NULL, x = "Depth (ft)")
dev.off()

