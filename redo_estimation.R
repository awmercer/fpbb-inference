library(BART)
library(tidyverse)
library(bayesboot)
library(bestimate)
library(timefactory)



get_estimate_posteriors = function(samp_id, samp, ref, synth_pop_ids, x_vars, y_vars, draws, cores) {
  
  sp_wts = tibble(ids = synth_pop_ids) %>% group_by(ids) %>% 
    summarise(wt = n()) %>%
    pull(wt)
  
  from_start = timefactory()
  
  x_ref = ref[, x_vars]
  x_samp = samp[, x_vars]
  n_samp = nrow(samp)


  # Get subsample of reference
  ref_subsamp_ids = synth_pop_ids[sample(seq_along(synth_pop_ids), n_samp, replace = FALSE)]
  ref_subsamp = ref[ref_subsamp_ids, ]
  x_ref_subsamp = ref_subsamp[ , x_vars]
  
  ## Estimate response propensities
  origin = c(rep(1, n_samp), rep(0, n_samp))
  comb = bind_rows(x_samp, x_ref_subsamp)

  # Pre-filled in BART call with standard parameters
  bart_partial = partial(pbart2, ndpost = draws, 
                        verbose = TRUE, 
                        keeptrainfits = FALSE,
                        mc.cores = cores)
  
  cat("Fitting propensities ")
  propensity_timer = timefactory()
  propensity_fit = bart_partial(x.train = comb, 
                          y.train = origin)
                        
  sample_propensities = pbart_posterior(propensity_fit, 
                                        newdata = x_samp, 
                                        mc.cores = cores)
  cat(sprintf("%.1f\n", propensity_timer()))
  

  # Get OR fits - confounded and unconfounded
  
  y_fits_timer = timefactory()
  cat("Fitting y models ")
  y_fits_confounded = y_vars %>% 
    map(~bart_partial(x.train = x_samp, y.train = samp[[.x]]))

  y_fits_unconfounded = y_vars %>%
    map(~bart_partial(x.train = x_ref_subsamp, y.train = ref_subsamp[[.x]]))
  
  cat(sprintf("%.1f\n", y_fits_timer()))
  

  # Add posterior mean propensity score to x_comb/x_ref
  x_samp_prop = x_samp %>%
    mutate(log_pi_hat = log(rowMeans(sample_propensities)))
  
  dr_fits_timer = timefactory()
  cat("Fitting DR models ")
  y_psc_fits = y_vars %>%
    map(~bart_partial(x.train = x_samp_prop, y.train = samp[[.x]]))
  cat(sprintf("%.1f\n", dr_fits_timer()))
  
  cat("Saving BART fits ")
  save_timer = timefactory()
  # Save BART fits to file for reuse later
  saveRDS(file = sprintf("data/bart_models/bart_fits_%s.RDS", samp_id),
          object = list(sample_id = samp_id, 
              propensity_fit = propensity_fit,
              y_fits_confounded = y_fits_confounded,
              y_fits_unconfounded = y_fits_unconfounded,
              y_psc_fits = y_psc_fits)
  )
  cat(sprintf("%.1f\n", save_timer()))
  
  ## Estimate posteriors and other quantities
  
  est_timer = timefactory()
  cat("Starting estimates:\n")
  
  # Calculate weights as odds of being in the population over sample
  sample_weights = map(sample_propensities, ~(1 - .x)/.x)
  
  # Bayesian bootstrap weights to accomodate sampling variance
  bb_weights = t(rudirichlet(draws, n_samp) * n_samp) %>% 
    as_tibble() %>%
    as.list() %>%
    set_names(sprintf("bb_wt_%s", seq_along(.)))
  
  ## cross bb and sample weights
  draw_bb = cross(list(draw = seq_len(draws), bb=seq_len(draws)))  %>% transpose()
  
  res = list()

  # Estimate propensity weighted means
  res$y_bar_propwt = y_vars %>% map_dfc(function(y_var) {
    pmap_dbl(draw_bb, function(draw, bb) {
      weighted.mean(samp[[y_var]], sample_weights[[draw]] * bb_weights[[bb]])
      })
  })
  
  cat(sprintf("finished propensity means %.1f\n", est_timer()))
  
  # Estimate simple unweighted bayes bootstrap means
  res$y_bar_samp_bayesboot = map_dfc(y_vars, function(y_var) {
    map_dbl(bb_weights, ~weighted.mean(samp[[y_var]], .x))
  })
  cat(sprintf("finished bayesboot means %.1f\n", est_timer()))
  
  # Estimate basic OR means
  res$y_bar_pred = map_dfc(y_fits_confounded, function(y_fit) {
    y_hat_pos = pbart_posterior(y_fit, newdata = x_ref, mc.cores = cores)
    
    map_dbl(y_hat_pos, ~weighted.mean(.x, sp_wts))
    
  })
  
  cat(sprintf("finished pred means %.1f\n", est_timer()))
  

  # Estimate DR-RBC means
  res$y_bar_drrbc = map_dfc(y_vars, function(y_var) {
    # Get OR model for y_var
    pred_fit = y_fits_confounded[[y_var]]
    
    # Get posterior predicted values for sample based on OR model
    y_hat_pos_samp = pbart_posterior(pred_fit, newdata = x_samp, mc.cores = cores)
    
    # Get the posterior distribution for the OR mean based on ref
    y_bar_pred_pos = res$y_bar_pred[[y_var]]
    
    # For each posterior draw calculate a weighted mean residual
    # for each of the bb_weights times the propensity weight for
    # the draw.
    # Add this weighted mean residual to the OR estimate
    # for that draw
    pmap_dbl(draw_bb, function(draw, bb) {
      y_bar_pred_draw = y_bar_pred_pos[[draw]]
      resid_draw = samp[[y_var]] - y_hat_pos_samp[[draw]]
      cmb_wt = sample_weights[[draw]] * bb_weights[[bb]]
      wtd_mean_resid = weighted.mean(resid_draw, cmb_wt)
      y_bar_drrbc = y_bar_pred_draw + wtd_mean_resid
    })
  })
  cat(sprintf("finished DR RBC means %.1f\n", est_timer()))
  
  # Esitimate OR-PSC means
  x_ref_prop = x_ref %>% 
    mutate(log_pi_hat = log(pbart_posterior(propensity_fit, newdata = x_ref, mc.cores = cores, return_posterior_mean = TRUE)))
  
  res$y_bar_drpsc = map_dfc(y_psc_fits, function(y_fit) {
    pos = pbart_posterior(y_fit, newdata = x_ref_prop, mc.cores = cores) 
    map_dbl(pos, ~weighted.mean(.x, sp_wts))
  })
  cat(sprintf("finished DR PSC means %.1f\n", est_timer()))
  cat(sprintf("Finished everything %.1f\n", from_start()))
  return(bind_rows(res, .id="est"))
}

np = readRDS("data/cleaned/cleaned_np_civic_data.RDS") 
cps = readRDS("data/cleaned/cps_civic_full_edited.RDS") 


x_vars = c("age", "sex", "racethn", "educcat", "fcregion")
y_vars = str_subset(names(np), "y_") %>% set_names()
np_samples = unique(np$sample_id)

synth_pop_ids = fpbb_synth_pops(weights = cps$pwsrwgt, 
                                L=1, 
                                N = nrow(cps) * 100, 
                                return_weights = FALSE)[[1]]

draws = 1000
est_pos_full = map(np_samples, function(samp_id) {
  est_pos = get_estimate_posteriors(samp_id = samp_id, 
                                    samp = filter(np, sample_id == samp_id), 
                                    ref = cps, 
                                    synth_pop_ids = synth_pop_ids, 
                                    x_vars = x_vars, 
                                    y_vars = y_vars, 
                                    draws = draws, 
                                    cores = 2)
  saveRDS(est_pos, sprintf("data/posteriors/est_pos_%s.RDS", samp_id))
})
saveRDS(est_pos_full, "data/posteriors/est_pos_full.RDS")

