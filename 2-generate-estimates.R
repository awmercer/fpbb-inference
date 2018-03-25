library(BART)
library(tidyverse)
library(bayesboot)
library(bestimate)
library(timefactory)
library(stringr)

## NOTE: timefactory and bestimate can be installed with:
## devtools::install_github("awmercer/timefactory")
## devtools::install_github("awmercer/bestimate")

# timefactory is for timing code
# bestimate contains functions to make working with BART easier

get_estimate_posteriors = function(samp_id,
                                   samp,
                                   ref,
                                   synth_pop_ids,
                                   x_vars,
                                   y_vars,
                                   draws,
                                   pweight_synth_pops,
                                   cores) {
  from_start = timefactory()
  
  
  ## Convert synthetic population idices into frequency weights for each
  ## record in the synthetic population
  sp_wts = tibble(ids = synth_pop_ids) %>% group_by(ids) %>%
    summarise(wt = n()) %>%
    pull(wt)
  
  
  # List containing output
  res = list()
  
  ## Convenience data structures
  x_ref = ref[, x_vars]
  x_samp = samp[, x_vars]
  n_samp = nrow(samp)
  
  
  # Get subsample of reference
  ref_subsamp_ids = synth_pop_ids[sample(seq_along(synth_pop_ids), n_samp, replace = FALSE)]
  ref_subsamp = ref[ref_subsamp_ids,]
  x_ref_subsamp = ref_subsamp[, x_vars]
  
  ## Estimate response propensities
  origin = c(rep(1, n_samp), rep(0, n_samp))
  comb = bind_rows(x_samp, x_ref_subsamp)
  
  # Pre-filled in BART call with standard parameters
  bart_partial = partial(
    pbart2,
    ndpost = draws,
    verbose = FALSE,
    keeptrainfits = FALSE,
    mc.cores = cores,
    nskip = 1000
  )
  
  cat("Fitting propensities ")
  propensity_timer = timefactory()
  propensity_fit = bart_partial(x.train = comb,
                                y.train = origin)
  
  sample_propensities = pbart_posterior(propensity_fit,
                                        newdata = x_samp,
                                        mc.cores = cores)
  cat(sprintf("%.1f\n", propensity_timer()))
  
  
  # Fit OR models - confounded and unconfounded
  
  y_fits_timer = timefactory()
  cat("Fitting y models ")
  y_fits_confounded = y_vars %>%
    map( ~ bart_partial(x.train = x_samp, y.train = samp[[.x]]))
  
  y_fits_unconfounded = y_vars %>%
    map( ~ bart_partial(x.train = x_ref_subsamp, y.train = ref_subsamp[[.x]]))
  
  cat(sprintf("%.1f\n", y_fits_timer()))
  
  
  # Add posterior mean propensity score to x_samp for OR-PSC
  x_samp_prop = x_samp %>%
    mutate(pi_hat = rowMeans(sample_propensities))
  
  dr_fits_timer = timefactory()
  cat("Fitting OR-PSC models ")
  y_psc_fits = y_vars %>%
    map( ~ bart_partial(x.train = x_samp_prop, y.train = samp[[.x]]))
  cat(sprintf("%.1f\n", dr_fits_timer()))
  
  cat("Saving BART fits ")
  save_timer = timefactory()
  # Save BART fits to file for reuse later
  saveRDS(
    file = sprintf("data/bart_models/bart_fits_%s.RDS", samp_id),
    object = list(
      sample_id = samp_id,
      propensity_fit = propensity_fit,
      y_fits_confounded = y_fits_confounded,
      y_fits_unconfounded = y_fits_unconfounded,
      y_psc_fits = y_psc_fits,
      synth_pop_ids = synth_pop_ids
    )
  )
  cat(sprintf("%.1f\n", save_timer()))
  
  ## Estimate posteriors and other quantities
  
  est_timer = timefactory()
  cat("Starting estimates:\n")
  
  # Calculate weights as odds of being in the population over sample
  sample_weights = map(sample_propensities, ~ (1 - .x) / .x)
  
  # For each propensity weight create a set of FPBB weights
  cat("Creating fpbb propensity weights ")
  fpbb_timer = timefactory()
  sink("/dev/null")
  sample_weight_synth_pops = map(
    sample_weights,
    ~ fpbb_synth_pops(
      weights = .x,
      L = pweight_synth_pops,
      N = length(.x) * 20
    )
  )
  sink()
  cat(sprintf("%.1f\n", fpbb_timer()))
  
  # Get "true" population means
  res$y_bar_pop = y_vars %>%
    map_dfc(function(y_var) {
      weighted.mean(ref[[y_var]], sp_wts)
    })

  # Estimate propensity weighted means
  res$y_bar_propwt = y_vars %>% map_dfc(function(y_var) {
    map(sample_weight_synth_pops, function(sp_wts) {
      map_dbl(sp_wts, function(wt) {
        weighted.mean(samp[[y_var]], wt)
      })
    }) %>% unlist()
  })
  
  cat(sprintf("finished propensity means %.1f\n", est_timer()))
  
  # Bayesian bootstrap weights to simulate SRS sampling variance
  bb_weights = t(rudirichlet(draws, n_samp) * n_samp) %>%
    as_tibble() %>%
    as.list() %>%
    set_names(sprintf("bb_wt_%s", seq_along(.)))
  
  # Estimate simple unweighted bayes bootstrap means
  res$y_bar_samp_bayesboot = map_dfc(y_vars, function(y_var) {
    map_dbl(bb_weights, ~ weighted.mean(samp[[y_var]], .x))
  })
  cat(sprintf("finished bayesboot means %.1f\n", est_timer()))
  
  # Estimate basic OR means
  res$y_bar_pred = map_dfc(y_fits_confounded, function(y_fit) {
    y_hat_pos = pbart_posterior(y_fit, newdata = x_ref, mc.cores = cores)
    
    map_dbl(y_hat_pos, ~ weighted.mean(.x, sp_wts))
    
  })
  
  cat(sprintf("finished pred means %.1f\n", est_timer()))

  # Estimate DR-RBC means
  res$y_bar_drrbc = map_dfc(y_vars, function(y_var) {
    # Get the posterior distribution for the OR mean based on ref
    y_bar_pred_pos = res$y_bar_pred[[y_var]]
    
    # Get OR model for y_var
    pred_fit = y_fits_confounded[[y_var]]
    
    # Get posterior predicted values for sample based on OR model
    y_hat_pos_samp = pbart_posterior(pred_fit, newdata = x_samp, mc.cores = cores)
    
    # Calculate the OR-RBC mean for each sp weight associted with each
    # posterior draw
    pmap(list(y_bar_pred_pos, y_hat_pos_samp, sample_weight_synth_pops),
         function(y_bar, y_hat, sp_wts) {
           resid = samp[[y_var]] - y_hat
           
           # For each sp_weight associated with the draw
           # calculate a weighted mean residual and add it
           # to the predicted mean for that draw
           map_dbl(sp_wts, function(wt) {
             y_bar + weighted.mean(resid, wt)
           }) %>% unlist()
           
         }) %>% unlist()
  })
  cat(sprintf("finished DR RBC means %.1f\n", est_timer()))
  
  # Get propensities for reference sample
  ref_propensities = pbart_posterior(propensity_fit, newdata = x_ref, mc.cores = cores)

  # Esitimate OR-PSC means
  x_ref_prop = x_ref %>%
    mutate(pi_hat = rowMeans(ref_propensities))
  
  res$y_bar_drpsc = map_dfc(y_psc_fits, function(y_fit) {
    pos = pbart_posterior(y_fit, newdata = x_ref_prop, mc.cores = cores)
    map_dbl(pos, ~ weighted.mean(.x, sp_wts))
  })
  
  cat(sprintf("finished DR PSC means %.1f\n", est_timer()))
  
  # Estimate quantities for bias decomposition
  ref_phi = map2(sample_propensities, ref_propensities,
                 function(s_prop, ref_prop) {
                   min_s_prop = min(s_prop)
                   phi = ref_prop >= min_s_prop
                 })

  res$y_bar_samp_confounded = y_vars %>%
    map_dfc(function(y_var) {
      y_pos = pbart_posterior(y_fits_confounded[[y_var]],
                              newdata = x_samp,
                              mc.cores = cores)
      y_bar_samp_confounded = colMeans(y_pos)
    })
  
  res$y_bar_samp_unconfounded = y_vars %>%
    map_dfc(function(y_var) {
      y_pos = pbart_posterior(y_fits_unconfounded[[y_var]],
                              newdata = x_samp,
                              mc.cores = cores)
      y_bar_samp_unconfounded = colMeans(y_pos)
    })
  
  # Unconfounded estimates for full population
  y_bar_pop = y_vars %>%
    map(function(y_var) {
      y_pos = pbart_posterior(y_fits_unconfounded[[y_var]],
                              newdata = x_ref,
                              mc.cores = cores)
      list(
        # Posterior for for unconfounded population mean
        y_bar_pop_unconfounded = colMeans(y_pos),
        
        # Posterior for unconfounded population mean among region
        # of common support
        y_bar_pop_unconfounded_cs = map2_dbl(y_pos, ref_phi, 
                                             function(y, phi) {
                                               weighted.mean(y, phi)
                                             })
      )
      
    }) %>%
    transpose()
  
  res = c(res, y_bar_pop)
    
  cat(sprintf("Finished everything %.1f\n", from_start()))
  return(bind_rows(res, .id = "est"))
}

np = readRDS("data/cleaned/cleaned_np_civic_data.RDS")
cps = readRDS("data/cleaned/cps_civic_full_edited.RDS")
draws = 1000
pweight_synth_pops = 25
save_output = TRUE

x_vars = c("age", "sex", "racethn", "educcat", "fcregion")
y_vars = str_subset(names(np), "y_") %>% set_names()
np_samples = unique(np$sample_id) %>% set_names()


## Comment out when not testing
# np = filter(np, sample_id %in% c("A", "B")) %>% sample_n(400)
# cps = sample_n(cps, 500)
# draws = 10
# pweight_synth_pops = 10
# save_output = FALSE
# np_samples = np_samples[1:2]
# y_vars = y_vars[1:2]
# save_output = FALSE
# np_samples = "A"
# y_vars = y_vars[1]


# Create synthetic population for use as reference sample
set.seed(1234)
synth_pop_ids = fpbb_synth_pops(
  weights = cps$pwsrwgt,
  L = 1,
  N = nrow(cps) * 100,
  return_weights = FALSE
)

## Loop over each sample and estimate all of the necessary conditional means
start_timer = timefactory()
est_pos_full = np_samples %>%
  set_names() %>%
  map(function(samp_id) {
    est_pos = get_estimate_posteriors(
      samp_id = samp_id,
      samp = filter(np, sample_id == samp_id),
      ref = cps,
      synth_pop_ids = synth_pop_ids[[1]],
      x_vars = x_vars,
      y_vars = y_vars,
      draws = draws,
      pweight_synth_pops = pweight_synth_pops,
      cores = 10
    )
    if (save_output) {
      saveRDS(est_pos,
              sprintf("data/posteriors/est_pos_%s.RDS", samp_id))
    }
    est_pos
  }) %>% bind_rows(.id = "sample_id")

if (save_output) {
  saveRDS(est_pos_full, "data/posteriors/est_pos_full.RDS")
}
cat(sprintf("Whole thing took %.1f seconds.\n", start_timer()))


### Get minimum inclusion propensitities for each sample and
### calculate the portion of the population with common support

synth_pop_wts = synth_pop_ids %>% group_by(sp_idx_1) %>%
  arrange(sp_idx_1) %>%
  summarise(sp_wt = n()) %>% pull(sp_wt)

fit_files = list.files("data/bart_models", full.names = TRUE)

pop_common_support = map(np_samples, function(samp_id) {
  samp = np %>% filter(sample_id == samp_id)
  
  fits = readRDS(sprintf("data/bart_models/bart_fits_%s.RDS", samp_id))
  
  # Get % phi for population based on propensity model
  samp_propensities = pbart_posterior(fits$propensity_fit,
                                      newdata = samp,
                                      mc.cores = 10)
  samp_mins = map_dbl(samp_propensities, min)
  pop_propensities = pbart_posterior(fits$propensity_fit,
                                     newdata = cps,
                                     mc.cores = 10)
  
  pop_phi = map2_dfc(pop_propensities, samp_mins, ~ .x < .y) %>%
    colMeans()
  
  tibble(pct_common_support = pop_phi,
         samp_min_pi = samp_mins)
  
}) %>% bind_rows(.id = "sample_id")

saveRDS(pop_common_support, "data/posteriors/common_support.RDS")
