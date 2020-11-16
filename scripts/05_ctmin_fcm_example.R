######################################################################
######################################################################
######################################################################
# - Thermal limit studies and sample size 
# - Script 005:  FCM diet case-study
# - Author:      Guy F. Sutton
# - Affiliation: Center for Biological Control, 
#                Rhodes University, South Africa
# - Contact:     g.sutton@ru.ac.za
# - Date:        09/11/2020
######################################################################
######################################################################
######################################################################

######################################################################
# Load packages ------------------------------------------------------
######################################################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, 
               tidyr,
               readr,
               retrodesign)

# Load ThermalSampleR package
# To download on first use, remove # from line below 
# devtools::install_github("CJMvS/ThermalSampleR")
library(ThermalSampleR)

######################################################################
# Set global defaults ------------------------------------------------
######################################################################

# Set ggplot theme (makes nice plots)
theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", 
                                              fill = NA),
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), 
                                                            "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), 
                                                            "mm")),
                  legend.position = "none"))

######################################################################
# Load raw bootstrap resamples ---------------------------------------
######################################################################

# We have already drawn our bootstrap samples, load that data 
# - This is the output from script 01. 

# Load raw data 
boot_tci <- readr::read_csv("./data_clean/ct_min_bootstrap_with_ci.csv")
head(boot_tci)

# Make insect_sp column into a factor 
boot_tci <- boot_tci %>%
  dplyr::mutate(insect_sp = as.factor(insect_sp))
head(boot_tci)

# Also, load raw data (without CI's)
raw_data <- readr::read_csv("./data_clean/ct_min_data_clean.csv")
head(raw_data)

# Make insect_sp column into a factor 
raw_data <- raw_data %>%
  dplyr::mutate(col = as.factor(insect_sp))
head(raw_data)

###########################################################################
# Use ThermalSampleR to evaluate sample size ------------------------------
###########################################################################

# Subset data 
raw_data_ex <- raw_data %>%
  dplyr::filter(insect_sp %in% c("Thaumatotibia leucotreta_std_diet_CO_30",
                                 "Thaumatotibia leucotreta_tre_diet_CO_30"))

# Let's assume we only had the first 20 samples (n = 20), for each species
raw_data_ex <- raw_data_ex %>%
  dplyr::filter(between(sample_size, 1, 20))
head(raw_data_ex)

# Perform a t-test to compare means 
t.test(ct_val ~ insect_sp, 
       data = raw_data_ex, 
       var.equal = FALSE)

# A significant p-value - this is what we were hoping for! 

# Let's use ThermalSampleR to see whether our estimate is valid
# and whether additional sampling is likely to yield a more precise estimate. 

# Perform bootstrap re-sampling
# - Set seed to make re-sampling reproducible
set.seed(2012)
sims_20 <- ThermalSampleR::boot_two(data = raw_data_ex,
                                    groups_col = col,
                                    group1 = "Thaumatotibia leucotreta_std_diet_CO_30",
                                    group2 = "Thaumatotibia leucotreta_tre_diet_CO_30",
                                    response= ct_val,
                                    n_max = 50,
                                    iter = 499)

# Plot bootstrap samples 
plots <- ThermalSampleR::plot_two_groups(x = sims_20,
                                         n_min = 3,
                                         n_max = 20,
                                         legend.position = c(0.7, 0.85)) 
plots

# Extract mean CTmin values reported in paper for n = 20 
fcm_1 <- raw_data_ex %>%
  dplyr::filter(col == "Thaumatotibia leucotreta_std_diet_CO_30") %>%
  dplyr::summarise(mean_ct = mean(ct_val),
                   n = n(),
                   sd_ct = sd(ct_val),
                   se_ct = sd_ct/sqrt(n))
fcm_1

fcm_2 <- raw_data_ex %>%
  dplyr::filter(col == "Thaumatotibia leucotreta_tre_diet_CO_30") %>%
  dplyr::summarise(mean_ct = mean(ct_val),
                   n = n(),
                   sd_ct = sd(ct_val),
                   se_ct = sd_ct/sqrt(n))
fcm_2

# Extract difference in CTmin values reported in paper for n = 20 
diff_20 <- sims_20 %>%
  dplyr::filter(sample_size == 20) %>%
  dplyr::select(sample_size, lower_ci, upper_ci, mean_diff) %>%
  dplyr::summarise(lower_ci = mean(lower_ci),
                   upper_ci = mean(upper_ci),
                   mean_diff = mean(mean_diff))
diff_20

# Second step, let's evaluate our study using design analysis calculations.
# - Following Gelman and Carlin (2014), calculate:
# (1) Type S errors - probability we got the sign of the relationship wrong.
# (2) Type M errors - factor by which we may overestimate difference in means. 

# We need pooled standard error 
raw_data_calc <- raw_data_ex %>%
  group_by(insect_sp) %>%
  mutate(ct_mean = mean(ct_val), 
         sd_mean = sd(ct_val),
         n = n()) %>%
  select(-sample_size, ct_val) %>%
  slice(1) %>%
  ungroup()
raw_data_calc

# Now, let's define our parameters 
# - n1 and n2 = samples sizes per grou
# - s21 and s2 = standard deviation of CTmin per group
n1 <- raw_data_calc %>%
  slice(1) %>%
  pull(n)
n2 <- raw_data_calc %>%
  slice(2) %>%
  pull(n)
s1 <- raw_data_calc %>%
  slice(1) %>%
  pull(sd_mean)
s2 <- raw_data_calc %>%
  slice(2) %>%
  pull(sd_mean)

# Finally, calculate pooled standard error
first_term <- (n1-1)*s1^2
second_term <- (n2-1)*s2^2
pooled_df <- (n1 + n2)-2
s_pooled <- sqrt((first_term + second_term)/pooled_df)
se_diff <- s_pooled * sqrt((1/n1)+(1/n2))  
se_diff

# Okay, now evaluate our study by Gelman and Carlin's 
# type S and type M errors by reading off the graph 

# Run simulations across a range of possible effect sizes and pooled sd
possible_effects <- c(0.03, 1.3)

# Plug in the se_diff value after the 'possible effects, x))
(effect_pairs_00 <- data.frame(possible_effects, 
                               type_s(possible_effects,se_diff),
                               type_m(possible_effects,se_diff)))

# If we just needed the FCM diet treatment CTmin's to be statistically different, 
# we could have been very confident, with a ~ 0% chance that we got the sign
# of the effect wrong. 
# However, we could be overestimating the different in means 
# by a factor of 42x.

# - The ThermalSampleR simulations show us that the width of the 
#   CI for difference in means is still declining with additional samples added,
#   and is only predicted to approach an asymptote around ~ n = 30. 

# In the next section, we will analyse the full dataset (n = 30),
# simulating an experimental design where we tested n = 20, performed the 
# analyses above (where we estimate that n = 30 may give us a more valid
# assessment of CTMin for our study and be statistically valid). 

# Subset data 
raw_data_ex <- raw_data %>%
  dplyr::filter(insect_sp %in% c("Thaumatotibia leucotreta_std_diet_CO_30",
                                 "Thaumatotibia leucotreta_tre_diet_CO_30"))

# Use all samples for each species (n = 30)
raw_data_ex <- raw_data_ex %>%
  dplyr::filter(between(sample_size, 1, 30))
head(raw_data_ex)

# Perform a t-test to compare means 
t.test(ct_val ~ insect_sp, 
       data = raw_data_ex, 
       var.equal = FALSE)

# A significant p-value - this is what you were hoping for! 

# Let's use ThermalSampleR to see whether our estimate is valid
# and whether additional sampling is likely to yield a more precise estimate. 

# Perform bootstrap resampling
sims_30 <- ThermalSampleR::boot_two(data = raw_data_ex,
                                    groups_col = col,
                                    group1 = "Thaumatotibia leucotreta_std_diet_CO_30",
                                    group2 = "Thaumatotibia leucotreta_tre_diet_CO_30",
                                    response= ct_val,
                                    n_max = 50,
                                    iter = 499)

# Plot bootstrap samples 
plots <- plot_two_groups(x = sims_30,
                         n_min = 3,
                         n_max = 30,
                         legend.position = c(0.7, 0.85)) 
plots

# Extract mean CTmin values reported in paper for n = 30 
fcm_1 <- raw_data_ex %>%
  dplyr::filter(col == "Thaumatotibia leucotreta_std_diet_CO_30") %>%
  dplyr::summarise(mean_ct = mean(ct_val),
                   n = n(),
                   sd_ct = sd(ct_val),
                   se_ct = sd_ct/sqrt(n))
fcm_1

fcm_2 <- raw_data_ex %>%
  dplyr::filter(col == "Thaumatotibia leucotreta_tre_diet_CO_30") %>%
  dplyr::summarise(mean_ct = mean(ct_val),
                   n = n(),
                   sd_ct = sd(ct_val),
                   se_ct = sd_ct/sqrt(n))
fcm_2

# Extract difference in CTmin values reported in paper for n = 30 
diff_30 <- sims_30 %>%
  dplyr::filter(sample_size == 30) %>%
  dplyr::select(sample_size, lower_ci, upper_ci, mean_diff) %>%
  dplyr::summarise(lower_ci = mean(lower_ci),
                   upper_ci = mean(upper_ci),
                   mean_diff = mean(mean_diff))
diff_30

# Second step, let's evaluate our study using design analysis calculations.
# - Following Gelman and Carlin (2014), calculate:
# (1) Type S errors - probability we got the sign of the relationship wrong.
# (2) Type M errors - factor by which we may overestimate difference in means. 

# We need pooled standard error 
raw_data_calc <- raw_data_ex %>%
  group_by(insect_sp) %>%
  mutate(ct_mean = mean(ct_val), 
         sd_mean = sd(ct_val),
         n = n()) %>%
  select(-sample_size, ct_val) %>%
  slice(1) %>%
  ungroup()
raw_data_calc

# Now, let's define our parameters 
n1 <- raw_data_calc %>%
  slice(1) %>%
  pull(n)
n2 <- raw_data_calc %>%
  slice(2) %>%
  pull(n)
s1 <- raw_data_calc %>%
  slice(1) %>%
  pull(sd_mean)
s2 <- raw_data_calc %>%
  slice(2) %>%
  pull(sd_mean)

first_term <- (n1-1)*s1^2
second_term <- (n2-1)*s2^2
pooled_df <- (n1 + n2)-2
s_pooled <- sqrt((first_term + second_term)/pooled_df)
se_diff <- s_pooled * sqrt((1/n1)+(1/n2))  

# Finally, let's see what our pooled se is
se_diff

# Okay, now evaluate our study by Gelman and Carlin's 
# type S and type M errors by reading off the graph 

# Run simulations across a range of possible effect sizes and s.e.
possible_effects <- c(0.15, 1.1)

# Plug in the se_diff value after the 'possible effects, x))
(effect_pairs_00 <- data.frame(possible_effects, 
                               type_s(possible_effects,se_diff),
                               type_m(possible_effects,se_diff)))

# If we just needed the FCM diet treatment CTmin's to be different, 
# we could have been very confident, with a ~0% chance that we got the sign
# of the effect wrong. 

# However, we could be still overestimating the different in means 
# by a factor of 7x.

# - The ThermalSampleR simulations show us that the width of the 
#   CI for difference in means may still decline with additional samples added,
#   above n = 30.
# - Depending on the degree of precision and logistics of testing more 
#   individuals, researchers may want to pursue adding more samples... 