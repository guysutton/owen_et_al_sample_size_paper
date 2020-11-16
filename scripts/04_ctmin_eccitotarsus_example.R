######################################################################
######################################################################
######################################################################
# - Thermal limit studies and sample size 
# - Script 004:  Eccritotarsus case-study
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

# Install our ThermalSampleR package from GitHub
# devtools::install_github("CJMvS/ThermalSampleR")
library(ThermalSampleR)

###########################################################################
# Set global defaults -----------------------------------------------------
###########################################################################

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

#####################################################################
# Load raw bootstrap resamples --------------------------------------
#####################################################################

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
  dplyr::filter(insect_sp %in% c("Eccritotarsus catarinensis_QTG",
                                 "Eccritotarsus eichhorniae_QTG"))

# Let's assume we only had the first 15 samples (n = 15), for each species
raw_data_ex <- raw_data_ex %>%
  dplyr::filter(between(sample_size, 1, 15))
head(raw_data_ex)

# Perform a t-test to compare means 
t.test(ct_val ~ insect_sp, 
       data = raw_data_ex, 
       var.equal = FALSE)

# A significant p-value - this is what you were hoping for! 
# Difference is estimated to be somewhere between 1.1 and 3.6 degrees. 

# Let's assume that we needed the new Eccritotarsus species 
# to have a CTmin > 1.2 degrees lower to warrant studying further. 
# - Our first assessment of n = 15 samples above showed the  
#   difference is estimated to be somewhere between 1.1 and 3.6 degrees.

# Let's use ThermalSampleR to see whether our estimate is valid
# and whether additional sampling is likely to yield a more precise estimate. 

# Perform bootstrap resampling
# - Set seed to make re-sampling reproducible
set.seed(2012)
sims_15 <- ThermalSampleR::boot_two(data = raw_data_ex,
                                    groups_col = col,
                                    group1 = "Eccritotarsus catarinensis_QTG",
                                    group2 = "Eccritotarsus eichhorniae_QTG",
                                    response= ct_val,
                                    n_max = 50,
                                    iter = 499)

# Plot bootstrap samples 
plots <- plot_two_groups(x = sims_15,
                         n_min = 3,
                         n_max = 15,
                         legend.position = c(0.8, 0.2)) 
plots

# Extract mean CTmin values reported in paper for n = 15 
e_cat_15 <- raw_data_ex %>%
  dplyr::filter(col == "Eccritotarsus catarinensis_QTG") %>%
  dplyr::summarise(mean_ct = mean(ct_val),
                   n = n(),
                   sd_ct = sd(ct_val),
                   se_ct = sd_ct/sqrt(n))
e_cat_15

e_eic_15 <- raw_data_ex %>%
  dplyr::filter(col == "Eccritotarsus eichhorniae_QTG") %>%
  dplyr::summarise(mean_ct = mean(ct_val),
                   n = n(),
                   sd_ct = sd(ct_val),
                   se_ct = sd_ct/sqrt(n))
e_eic_15

# Extract difference in CTmin values reported in paper for n = 15 
diff_15 <- sims_15 %>%
  dplyr::filter(sample_size == 15) %>%
  dplyr::select(sample_size, lower_ci, upper_ci, mean_diff) %>%
  dplyr::summarise(lower_ci = mean(lower_ci),
                   upper_ci = mean(upper_ci),
                   mean_diff = mean(mean_diff))
diff_15

# Second step, let's evaluate our study using design analysis calculations.
# - Following Gelman and Carlin (2014), calculate:
# (1) Type S errors - probability we got the sign of the relationship wrong.
# (2) Type M errors - factor by which we may overestimate difference in means. 

# We need pooled standard error 
raw_data_calc <- raw_data_ex %>%
  dplyr::group_by(insect_sp) %>%
  dplyr::mutate(ct_mean = mean(ct_val), 
                sd_mean = sd(ct_val),
                n = n()) %>%
  dplyr::select(-sample_size, ct_val) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
raw_data_calc

# Now, let's define our parameters 
n1 <- raw_data_calc %>%
  dplyr::slice(1) %>%
  dplyr::pull(n)
n2 <- raw_data_calc %>%
  dplyr::slice(2) %>%
  dplyr::pull(n)
s1 <- raw_data_calc %>%
  dplyr::slice(1) %>%
  dplyr::pull(sd_mean)
s2 <- raw_data_calc %>%
  dplyr::slice(2) %>%
  dplyr::pull(sd_mean)

# Calculate pooled sd 
first_term <- (n1-1)*s1^2
second_term <- (n2-1)*s2^2
pooled_df <- (n1 + n2)-2
s_pooled <- sqrt((first_term + second_term)/pooled_df)
se_diff <- s_pooled * sqrt((1/n1)+(1/n2))  
se_diff

# Okay, now evaluate our study by Gelman and Carlin's 
# type S and type M errors by reading off the graph 

# Run simulations across a range of possible effect sizes and s.e.
possible_effects <- c(1.1, 3.6)

# Plug in the se_diff value after the 'possible effects, x))
(effect_pairs_00 <- data.frame(possible_effects, 
                              type_s(possible_effects,se_diff),
                              type_m(possible_effects,se_diff)))

# If we just needed the two Eccritotarsus CTmin's to be different, 
# we could have been very confident, with a ~ 0% chance that we got the sign
# of the effect wrong. 

# Inspecting the results shows us that the CTmin of E. eichhorniae is 
# actually HIGHER than that of E. catarinensis, so we have statistically
# robust findings at n = 15 individuals tested, that we don't need
# more testing to answer our study question. 
