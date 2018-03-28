library(tidyverse)

bart_list = list.files("data/bart_models", full.names = TRUE)


propensities = map(bart_list, function(bartfile) {
  fit = readRDS(bartfile)
  props = fit$propensity_fit$prob.train.mean
  n = length(props)/2
  
  tibble(sample_id = fit$sample_id,
         source = c(rep(fit$sample_id, n), 
                     rep("CPS", n)),
         pscore = props
  )
}) %>%
  bind_rows()

saveRDS(propensities, "data/propensities.RDS")
