results <- readRDS("~/Desktop/SUPR/results.rds")

methods <- c("NSP", "NOT", "MOSUM BUM", "MOSUM LP", "PELT", 
             "HSMUCE (0.1)", "HSMUCE (0.5)", 
             "MICH Ora Rev", "MICH AutoFast Rev")

uq_methods <- c("NSP", "MOSUM BUM", "MOSUM LP", 
                "HSMUCE (0.1)", "HSMUCE (0.5)", 
                "MICH Ora Rev", "MICH AutoFast Rev")

bias_star = results %>% 
  filter(method %in% methods) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
  group_by(T, L, min_space, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L, min_space) %>% 
  summarize(bias_star = min(bias),
            hausdorff_star = min(hausdorff, na.rm = TRUE),
            fpsle_star = min(fpsle, na.rm = TRUE),
            fnsle_star = min(fnsle, na.rm = TRUE))

bias_dag = results %>% 
  filter(method %in% methods & !grepl("Ora", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
  group_by(T, L, min_space, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L, min_space) %>% 
  summarize(bias_dag = min(bias),
            hausdorff_dag = min(hausdorff, na.rm = TRUE),
            fpsle_dag = min(fpsle, na.rm = TRUE),
            fnsle_dag = min(fnsle, na.rm = TRUE))

bias_ddag = results %>% 
  filter(method %in% uq_methods & !grepl("Ora", method)) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
  group_by(T, L, min_space, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
            hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE)) %>% 
  group_by(T, L, min_space) %>% 
  summarize(bias_ddag = min(bias),
            hausdorff_ddag = min(hausdorff, na.rm = TRUE),
            fpsle_ddag = min(fpsle, na.rm = TRUE),
            fnsle_ddag = min(fnsle, na.rm = TRUE))

results %>% 
  filter(method %in% methods) %>% 
  filter((delta == 0.5 | is.na(delta)) & (prior == 1e-3 | is.na(prior)) & (tol == 1e-7 | is.na(delta))) %>% 
  group_by(T, L, min_space, method) %>% 
  summarize(bias = mean(abs(L - L_est)),
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
  arrange(T, L, min_space, method) %>% 
  mutate(method = gsub("[a-z] ", "", method)) %>%
  mutate(bias = gsub(" ", "", paste0(bias_note, bias)),
         hausdorff = gsub(" ", "", paste0(haus_note, hausdorff)),
         fpsle = gsub(" ", "", paste0(fpsle_note, fpsle)), 
         fnsle = gsub(" ", "", paste0(fnsle_note, fnsle))) %>% 
  mutate_all(trimws) %>% 
  select(T, L, min_space, method, bias, hausdorff, fpsle, fnsle, ci_length, coverage, time) %>%
  write.csv("~/Desktop/SUPR/simulations/main_sim.csv", row.names = FALSE)
