######################################################################
######################################################################
######################################################################
# - Thermal limit studies and sample size 
# - Script 002:  One-sample simulations - confidence interval coverage
# - Author:      Guy F. Sutton
# - Affiliation: Centre for Biological Control, 
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
               data.table,
               viridis,
               Hmisc,
               broom)

###########################################################################
# Set global defaults -----------------------------------------------------
###########################################################################

# Set ggplot theme (makes nice plots)
theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill = NA),
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
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

#####################################################################
# Calculate coverage of pop. parameter (i.e. CTmin) -----------------
#####################################################################

# We have already calculated median CTmin values per species 
# (median for n = 30; max. sample size in our study) in script 01.
# Following Pearson and Groves (2013), we assumed that our maximum sample size (n = 30)
# is a reasonable approximation of the actual critical thermal limit value (population parameter)
# i.e. if we hypothetically had sampled the entire population
# Obviously, we must interpret these data with caution, as there is no possible way to 
# 100% accurately determine the population parameter, so we must derive an ESTIMATE.

# Coverage refers to the proportion of bootstrap resamples for which the 
# estimated population parameter falls within the bounds of a 95% confidence interval. 
# - Coverage provides an estimate of parameter accuracy. 

# Calculate the proportion of times the median_pop_val falls 
# within the lower_ci and upper_ci bounds (95%  CI)
bootstrap_cover <- boot_tci %>%
  dplyr::group_by(insect_sp, sample_size, iter) %>%
  dplyr::mutate(ci_falls = dplyr::case_when(
    median_pop_val < lower_ci ~ 0,
    median_pop_val > upper_ci ~ 0, 
    median_pop_val > lower_ci & median_pop_val < upper_ci ~ 1, 
    median_pop_val < upper_ci & median_pop_val > lower_ci ~ 1
  )) %>%
  dplyr::group_by(insect_sp, sample_size) %>%
  dplyr::mutate(n_correct = sum(ci_falls),
                max_n = max(iter),
                prop_correct = n_correct/max_n) %>%
  dplyr::slice(1)
head(bootstrap_cover)

#####################################################################
# Plot coverage of pop. parameter (i.e. CTmin), individual species --
#####################################################################

# Plot across all species 
ggplot(data = bootstrap_cover, aes(x = sample_size,
                                   y = prop_correct,
                                   colour = insect_sp)) + 
  geom_line() +
  scale_colour_grey() + 
  #scale_color_viridis(discrete = TRUE, option = "C") +
  labs(x = "Sample size (n)",
       y = "95% CI coverage\n (Proportion containing population parameter)",
       subtitle = " ") +
  scale_y_continuous(breaks = c(0.7, 0.75, 
                                0.8, 0.85, 0.9, 0.95, 1.00), 
                     limits = c(0.7, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dashed") + 
  facet_wrap(~ insect_sp)

# Save figure to file 
ggsave("./figures/figure_s1.png", 
       width = 20, 
       height = 12)

#####################################################################
# Plot coverage of pop. parameter (i.e. CTmin), pooled species --
#####################################################################

# Might be more clear by calculating coverage across all species 
# Manually calculate binomial confidence intervals
bin_prop_ci <- as_tibble(cbind(bootstrap_cover, 
                               Hmisc::binconf(bootstrap_cover$n_correct, 
                                              bootstrap_cover$max_n)))

# Extract binconf output 
aa <- as_tibble(bin_prop_ci$...15)

# Add binconf output back to bootstrap samples 
bin_prop_ci <- dplyr::bind_cols(bin_prop_ci, aa) %>%
  dplyr::select(insect_sp, 
                sample_size, 
                prop_correct = PointEst,
                lower_ci = Lower,
                upper_ci = Upper)

# Calculate summarised coverage vals and CI's across species 
bin_prop_ci %>%
  dplyr::group_by(sample_size) %>%
  dplyr::summarise(prop_correct = mean(prop_correct),
                   lower_ci = mean(lower_ci),
                   upper_ci = mean(upper_ci)) %>%
  ggplot(data = ., aes(x = sample_size,
                       y = prop_correct)) + 
  geom_ribbon(aes(ymin = lower_ci,
                  ymax = upper_ci),
              fill = "gray80") + 
  geom_line() +
  labs(x = "Sample size (n)",
       y = "95% CI coverage\n (Proportion containing population parameter)",
       subtitle = " ") +
  scale_y_continuous(breaks = c(0.8, 0.85, 0.9, 0.95, 1.00), 
                     limits = c(0.8, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_hline(yintercept = 0.90, linetype = "dotted")

# Save figure to file 
ggsave("./figures/ci_ctmin_ci_coverage_summary.png", 
       width = 6, 
       height = 6)


