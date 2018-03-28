library(tidyverse)
library(bestimate)
library(stringr)

cps = readRDS("data/cleaned/cps_civic_full_edited.RDS") %>%
  mutate(fcregion = factor(fcregion, levels = 1:4, labels=c("Northeast", "Midwest", "South", "West")))
np = readRDS("data/cleaned/cleaned_np_civic_data.RDS") %>%
  mutate(fcregion = factor(fcregion, levels = 1:4, labels=c("Northeast", "Midwest", "South", "West")))


set.seed(1234)
synth_pop_wts = fpbb_synth_pops(
  weights = cps$pwsrwgt,
  L = 1,
  N = nrow(cps) * 100,
  return_weights = TRUE
)
cps$sp_wt = synth_pop_wts$sp_wt_1

cps_y_means = cps %>%
  summarise_at(vars(starts_with("y")), weighted.mean, w = .$sp_wt) %>%
  gather(variable, CPS)

  
y_var_summary = np %>% 
  select(sample_id, starts_with("y")) %>%
  group_by(sample_id) %>%
  summarise_at(vars(starts_with("y")), mean) %>%
  gather(variable, value, -sample_id) %>%
  left_join(cps_y_means) %>%
  mutate(diff_from_cps = value - CPS) %>%
  select(-value) %>%
  spread(sample_id, diff_from_cps) %>%
  mutate_if(is.numeric, ~round(.x * 100)) %>%
  mutate_at(vars(-variable, -CPS), ~sprintf("%+.0f", .x)) %>%
  mutate(CPS = as.character(CPS)) %>%
  mutate(Variable = factor(variable, 
                           levels = c("y_always_vote_local",
                                      "y_talk_neighbor_weekly",
                                      "y_trust_neighbors_all_most",
                                      "y_civic_assoc_yes",
                                      "y_school_group_yes",
                                      "y_recreational_assoc_yes"),
                           labels = c("Always votes in local elections",
                                      "Talks to neighbors daily/weekly",
                                      "Trusts all/most people in neighborhood",
                                      "Participated in civic association",
                                      "Participated in school group",
                                      "Participated in recreational associaton")
                           )) %>%
  select(Variable, CPS, everything(), -variable) %>%
  mutate_all(~str_replace(.x, "[+-]0", "0"))


cps_demos = cps %>% mutate(Age = cut(age, 
                             breaks = c(-Inf, 24, 34, 44, 54, 64, Inf), 
                             labels = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+"))) %>%
  select(Age, Sex=sex, `Race/Ethnicity` = racethn, Education = `educcat`, Region=fcregion, sp_wt) %>%
  gather(Variable, Category, -sp_wt) %>%
  group_by(Variable, Category) %>%
  summarise(cat_sum = sum(sp_wt)) %>%
  group_by(Variable) %>%
  mutate(CPS = cat_sum/sum(cat_sum)) %>%
  select(-cat_sum)
                  

demo_summary = np %>% mutate(Age = cut(age, 
                                     breaks = c(-Inf, 24, 34, 44, 54, 64, Inf), 
                                     labels = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+"))) %>%
  select(sample_id, Age, Sex=sex, `Race/Ethnicity` = racethn, Education = `educcat`, Region=fcregion) %>%
  gather(Variable, Category, -sample_id) %>%
  group_by(sample_id, Variable, Category) %>%
  summarise(n = n()) %>%
  group_by(sample_id, Variable) %>%
  mutate(pct = n/sum(n)) %>%
  select(-n) %>%
  left_join(cps_demos) %>%
  mutate(diff_from_cps = pct - CPS) %>%
  select(-pct) %>%
  mutate_at(vars(CPS, diff_from_cps), ~round(.x * 100)) %>%
  mutate(CPS = as.character(CPS),
         diff_from_cps = sprintf("%+.0f", diff_from_cps),
         diff_from_cps = str_replace(diff_from_cps, "[+-]0", "0")) %>%
  spread(sample_id, diff_from_cps) %>%
  mutate(Category = factor(Category,
         levels = c("Male",
                    "Female",
                    "18-24",
                    "25-34",
                    "35-44",
                    "45-54",
                    "55-64",
                    "65+",
                    "HS or Less",
                    "Some College",
                    "College Grad",
                    "Non-Hispanic White",
                    "Non-Hispanic Black",
                    "Hispanic",
                    "Other",
                    "Northeast", 
                    "Midwest",
                    "South", 
                    "West")
                    
                    )) %>%
  arrange(Category)

saveRDS(y_var_summary, "data/y_var_summary.RDS")
saveRDS(demo_summary, "data/demo_summary.RDS")
                    

