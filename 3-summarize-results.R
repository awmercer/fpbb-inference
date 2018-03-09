library(tidyverse)
library(stringr)
library(glue)
library(forcats)

# lin_pos = list.files("data/posteriors", pattern = "lin_pos_[A-Z]", full.names = TRUE) %>%
#   set_names(str_extract(., "(?<=pos_).*?(?=\\.)")) %>%
#   map(readRDS) %>%
#   map(bind_rows) %>%
#   bind_rows() %>%
#   mutate(x_vars = "demos") %>%
#   rename(y_bar_wlsdr = ybar_wlsdr, y_bar_ols = ybar_ols)

full_pos = list.files("data/posteriors", pattern = "full_pos_[A-Z]", full.names = TRUE) %>%
  set_names(str_extract(., "(?<=pos_).*?(?=\\.)")) %>%
  map(readRDS) %>%
  bind_rows(.id = "sample_id") %>%
  mutate(x_vars = "demos") %>%
  filter(sp == "sp_wt_1")
# %>%
#   right_join(lin_pos) # full_pos had extra draws due to multicore

# See if multiple synthetic populations are necessary

estimates_long = full_pos %>%
  select(sample_id, y_var, sp, true_y_bar = y_bar_obs_ref_tot, 
         y_bar_propwt, 
         y_bar_pred, 
         y_bar_drwt, 
         y_bar_dr, 
         # y_bar_wlsdr, 
         # y_bar_ols, 
         y_bar_samp_bayesboot) %>%
  gather(est, value, starts_with("y_bar")) 

per_sp_posterior = estimates_long %>% group_by(sample_id, y_var, sp, est) %>%
  summarise(mean = mean(value),
            var = var(value),
            bias = mean(value - true_y_bar),
            rmse = sqrt(mean((value - true_y_bar)^2)))
            
all_sp_posterior = estimates_long %>% group_by(sample_id, y_var, est) %>%
  summarise(mean = mean(value),
            var = var(value),
            bias = mean(value - true_y_bar),
            rmse = sqrt(mean((value - true_y_bar)^2)))

variance_comparison = per_sp_posterior %>%
  group_by(sample_id, y_var, est) %>%
  summarise(avg_within_sp_var = mean(var)) %>%
  left_join(all_sp_posterior %>% select(sample_id, y_var, est, total_var = var)) %>%
  mutate(deff = total_var/avg_within_sp_var)

variance_comparison %>% group_by(est) %>%
  summarise(mean(deff))

  

# reg_pos = list.files("data/regularization_tests/posteriors", pattern = "pos_[A-Z]", full.names = TRUE) %>%
#   set_names(str_extract(., "(?<=pos_).*?(?=\\.)")) %>%
#   map(readRDS) %>%
#   bind_rows(.id = "sample_id") %>%
#   mutate(x_vars = "demos") 
# rename(full_pos, y_bar_ref_tot = tot_y_bar_ref)

summarize_pos = function(df) {
  
  df %>%
    rename(tot_y_bar_obs_ref = y_bar_obs_ref_tot) %>%
  gather(key = "est", value = "value", starts_with("y_bar")) %>%
  group_by(sample_id, x_vars, est, y_var) %>%
  mutate(y_bar_ref = mean(tot_y_bar_obs_ref),
         delta = value - y_bar_ref,
         ucl95 = quantile(value, 0.975),
         lcl95 = quantile(value, 0.025),
         covered = y_bar_ref >= lcl95 & y_bar_ref <= ucl95) %>%
  summarise(y_bar_ref = mean(y_bar_ref),
            pos_mean = mean(value),
            pos_sd = sd(value),
            bias = mean(delta),
            rmse = sqrt(mean(delta^2)),
            ucl95 = mean(ucl95),
            lcl95 = mean(lcl95),
            coverage = mean(covered)) %>%
  select(y_var, est, sample_id, everything())  
}





main_summary = full_pos %>% 
  select(sample_id, x_vars, y_var, sp, draw, y_bar_obs_ref_tot, 
         y_bar_samp_bayesboot, 
         y_bar_drwt, 
         y_bar_propwt, 
         y_bar_pred, 
         y_bar_dr, 
         y_bar_samp_bayesboot
         # y_bar_wlsdr,
         # y_bar_ols
         ) %>%
  summarize_pos()


main_summary %>% 
  filter(x_vars == "demos") %>%
  #filter(est != "y_bar_samp_bayesboot") %>%
  #filter(est %in% c("y_bar_drwt", "y_bar_propwt")) %>%
  ggplot(aes(y=abs(bias), x=sample_id, color = est)) +
  facet_wrap(~y_var) +
  geom_point(alpha = .7, size = 2, shape = 1, stroke = 1) +
  coord_flip() +
  theme_bw() +
  ggtitle("Bias summary")

main_summary %>% 
  filter(x_vars == "demos") %>%
  mutate(y_var = fct_reorder(y_var, rmse)) %>%
  #filter(est != "y_bar_samp_bayesboot") %>%
  #filter(est %in% c("y_bar_drwt", "y_bar_propwt")) %>%
  ggplot(aes(y=rmse, x=y_var, color = est)) +
  facet_wrap(~sample_id) +
  geom_point(alpha = .7) +
  coord_flip() +
  theme_bw() +
  ggtitle("RMSE summary")


by_sample = main_summary %>%
  group_by(sample_id, est) %>%
  summarise(mean_abs_bias = mean(abs(bias)),
            mean_rmse = mean(rmse),
            mean_var = mean(pos_sd)) 


by_variable = main_summary %>%
  group_by(y_var, est) %>%
  summarise(mean_abs_bias = mean(abs(bias)),
            mean_rmse = mean(rmse),
            mean_var = mean(pos_sd)) 


by_sample %>%
  mutate(est = fct_reorder(est, mean_abs_bias)) %>%
  ggplot(aes(y = mean_abs_bias, x=sample_id, color = est)) +
  geom_point(size = 3, shape=21, fill = "white", stroke = 1.5) +  
  scale_color_brewer(type = "qual", palette = 1) +
  theme_bw() +
  coord_flip() +
  ggtitle("Mean absolute bias by sample")

by_sample %>% 
  filter(est != "y_bar_samp_bayesboot") %>%
  arrange(sample_id, mean_var) %>% knitr::kable(digits = 3)

by_variable %>% 
  filter(est != "y_bar_samp_bayesboot") %>%
  arrange(mean_var) %>% slice(3:4)

main_summary %>%
  group_by(sample_id, est) %>%
  summarise(mean_rmse = mean(rmse)) %>%
  mutate(est = fct_reorder(est, mean_rmse)) %>%
  ggplot(aes(y = mean_rmse, x=sample_id, color = est)) +
  geom_point(size = 3, shape=21, fill = "white", stroke = 1.5) +  
  theme_bw() +
  coord_flip() +
  ggtitle("Mean RMSE by sample")


main_summary %>%
  group_by(y_var, est) %>%
  summarise(mean_abs_bias = mean(abs(bias))) %>%
  mutate(est = fct_reorder(est, mean_abs_bias)) %>%
  ggplot(aes(y = mean_abs_bias, x=y_var, color = est)) +
  geom_point(size = 3, shape=21, fill = "white", stroke = 1.5) +  
  theme_bw() +
  coord_flip() +
  ggtitle("Mean absolute bias by question")

main_summary %>%
  group_by(y_var, est) %>%
  summarise(mean_rmse = mean(rmse)) %>%
  mutate(est = fct_reorder(est, mean_rmse)) %>%
  ggplot(aes(y = mean_rmse, x=y_var, color = est)) +
  geom_point(size = 3, shape=21, fill = "white", stroke = 1.5) +  
  theme_bw() +
  coord_flip() +
  ggtitle("Mean RMSE by question")




pred_check = full_pos %>% 
  select(sample_id, y_var, matches("model_bias"), matches("model_rmse"), mean_drwt_term2) %>%
  select(-matches("_n?cs")) %>%
  group_by(sample_id, y_var) %>%
  summarise_if(is.numeric, mean)

pred_check %>% gather(est, value, -sample_id, -y_var) %>%
  filter(!str_detect(est, "rmse")) %>%
  ggplot(aes(y=value, x=sample_id, color=est)) +
  geom_point() +
  facet_wrap(~y_var) +
  coord_flip()

