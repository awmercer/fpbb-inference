library(tidyverse)
library(bestimate)
library(stringr)

np = readRDS("data/cleaned/cleaned_np_civic_data.RDS") 
cps = readRDS("data/cleaned/cps_civic_full_edited.RDS") 
sp_weights = readRDS("data/cleaned/sp_weights.RDS")

x_vars = c("age", "sex", "racethn", "educcat", "fcregion")
y_vars = str_subset(names(np), "y_")
np_samples = unique(np$sample_id)



### Posterior estimates for each sample
set.seed(1234)
walk(np_samples, function(np_id) {
  cat(sprintf("-----------------------Starting Sample %s----------------------\n", np_id))
  np_estimates = bestimate(samp = filter(np, sample_id==np_id),
                 ref = cps,
                 y_var_names = y_vars,
                 x_var_names = x_vars,
                 sp_wts = sp_weights,
                 propensity = TRUE,
                 prediction = TRUE,
                 double_robust_reg = TRUE,
                 double_robust_wt = TRUE,
                 dr_propensity_transform = log,
                 quality_measures = TRUE,
                 propensity_replicates = 1,
                 posterior_draws = 1000, 
                 bart_params = list(mc.cores = 13, ntree = 50))
  saveRDS(np_estimates, sprintf("data/posteriors/full_pos_%s.RDS", np_id))
})


