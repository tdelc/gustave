library(tidyverse)
library(gustave)

lfs_samp_area <- gustave::lfs_samp_area
lfs_samp_dwel <- gustave::lfs_samp_dwel 
lfs_samp_ind <- gustave::lfs_samp_ind

head(gustave::lfs_samp_dwel)

# Definition of a variance estimation wrapper
precision_ict <- qvar(
  
  # As before
  data = ict_sample,
  dissemination_dummy = "dissemination",
  dissemination_weight = "w_calib",
  id = "firm_id",
  scope_dummy = "scope",
  sampling_weight = "w_sample", 
  strata = "strata",
  nrc_weight = "w_nrc", 
  response_dummy = "resp", 
  hrg = "hrg",
  calibration_weight = "w_calib",
  calibration_var = c(paste0("N_", 58:63), paste0("turnover_", 58:63)),
  
  # Replacing the variables of interest by define = TRUE
  define = TRUE
  
)

precision_ict(ict_sample, mean(employees))


# CASE 1

lfs_samp_dwel$diss    <- TRUE
lfs_samp_dwel$w_dwel  <- 1 / lfs_samp_dwel$pik_dwel
lfs_samp_dwel$w_final <- 1 / lfs_samp_dwel$pik

qvar1 <- qvar_ms(
  data = lfs_samp_dwel,
  id = "id_dwel",
  dissemination_dummy = "diss",
  dissemination_weight = "w_final",
  sampling_weight = "w_final",
  define = TRUE
)

qvar1(lfs_samp_dwel, total(income))

qvar1b <- qvar(
  data = lfs_samp_dwel,
  id = "id_dwel",
  dissemination_dummy = "diss",
  dissemination_weight = "w_final",
  sampling_weight = "w_final",
  define = TRUE
)

qvar1b(lfs_samp_dwel, total(income))

# CASE 2

lfs_samp_area$w_area  <- 1 / lfs_samp_area$pik_area
lfs_samp_area$w_area  <- base::mean(lfs_samp_area$w_area)

lfs_samp_dwel <- lfs_samp_dwel %>%
  left_join(lfs_samp_area %>% select(id_area, w_area), by = "id_area")

lfs_samp_dwel$w_dwel  <- 1 / lfs_samp_dwel$pik_dwel
lfs_samp_dwel$w_final <- lfs_samp_dwel$w_dwel * lfs_samp_dwel$w_area

head(lfs_samp_area)
head(lfs_samp_dwel)

qvar2 <- qvar_ms(
  # One file per stage, from the highest to the lowest
  data = list(lfs_samp_area, lfs_samp_dwel),
  sampling_stages = 2,
  id = c("id_area", "id_dwel"),
  parent_id = c("id_area"),
  dissemination_dummy = "diss",
  dissemination_weight = "w_final",
  # Conditional sampling weights, one per stage
  sampling_weight = c("w_area", "w_dwel"),
  define = TRUE
)

qvar2(lfs_samp_dwel, total(income))


# CASE 3

lfs_samp_dwel$prep   <- 0.8
lfs_samp_dwel$w_nrc <- lfs_samp_dwel$w_dwel / lfs_samp_dwel$prep   # conditional!
lfs_samp_dwel$w_final <- lfs_samp_dwel$w_nrc * lfs_samp_dwel$w_area

head(lfs_samp_dwel)

qvar3 <- qvar_ms(
  data = list(lfs_samp_area, lfs_samp_dwel),
  sampling_stages = 2,
  id = c("id_area", "id_dwel"),
  parent_id = c("id_area"),
  dissemination_dummy = "diss",
  dissemination_weight = "w_final",
  sampling_weight = c("w_area", "w_dwel"),
  # Single names: assigned to the last stage by default
  nrc_weight = "w_nrc",
  # response_prob = "prep", # instead of nrc_weight
  response_dummy = "diss",
  define = TRUE
)
qvar3(lfs_samp_dwel, total(income))




# CASE 4

lfs_samp_dwel <- lfs_samp_dwel |> mutate(q_income = cut(income,5))
head(lfs_samp_dwel)

qvar4 <- qvar_ms(
  data = list(lfs_samp_area, lfs_samp_dwel),
  sampling_stages = 2,
  id = c("id_area", "id_dwel"),
  parent_id = c("id_area"),
  dissemination_dummy = "diss",
  dissemination_weight = "w_final",
  sampling_weight = c("w_area", "w_dwel"),
  # Single names: assigned to the last stage by default
  nrc_weight = "w_nrc",
  calibration_weight = "w_final",
  calibration_var = "q_income",
  response_dummy = "diss",
  define = TRUE
)
qvar4(lfs_samp_dwel, total(income))


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


