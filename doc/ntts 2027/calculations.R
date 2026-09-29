library(tidyverse)
library(gustave)

lfs_samp_area <- gustave::lfs_samp_area |> 
  mutate(
    strata = ifelse(id_area %in% c("A016","A007"), "A", "B"),
    w_area = 1 / pik_area,
    w_area = base::mean(w_area))

set.seed(123)

# Nombre de personnes au chomage par ménage
unemp_dwel <- gustave::lfs_samp_ind |> 
  group_by(id_dwel) |> 
  summarise(unemp = sum(unemp))

lfs_samp_dwel <- gustave::lfs_samp_dwel |> 
  left_join(lfs_samp_area %>% select(id_area, w_area), by = "id_area") |>  
  left_join(unemp_dwel, by = "id_dwel") |>  
  mutate(w_dwel = 1 / pik_dwel,
         w_final = w_dwel * w_area) |> 
  mutate(rand     = runif(n()),
         prep     = 0.8,
         rep      = rand <= prep,
         w_nrc    = w_dwel / prep,   # conditional!
         w_final  = ifelse(rep,w_nrc * w_area,NA),
         q_income = cut(income,5)
  )

wrapper_lfs <- qvar_ms(
  data = list(lfs_samp_area, lfs_samp_dwel),
  sampling_stages = 2,
  id = c("id_area", "id_dwel"),
  parent_id = c("id_area"),
  strata = list("strata",NULL),
  dissemination_dummy = "rep",
  dissemination_weight = "w_final",
  sampling_weight = c("w_area", "w_dwel"),
  nrc_weight = "w_nrc",
  response_dummy = "rep",
  calibration_weight = "w_final",
  calibration_var = "q_income",
  define = TRUE
)

inspect_wrapper(wrapper_lfs)

wrapper_lfs(lfs_samp_dwel, sum(unemp), diagnostics = T)

lfs_samp_dwel |> 
  ggplot() +
  aes(x = income) + 
  geom_histogram() + 
  facet_wrap(~id_area)


