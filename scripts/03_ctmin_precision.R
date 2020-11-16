######################################################################
######################################################################
######################################################################
# - Thermal limit studies and sample size 
# - Script 003:  One-sample simulations - confidence interval precision
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

# Load raw data (with CI's)
boot_tci <- readr::read_csv("./data_clean/ct_min_bootstrap_with_ci.csv")
head(boot_tci)

# Make insect_sp column into a factor 
boot_tci <- boot_tci %>%
  dplyr::mutate(insect_sp = as.factor(insect_sp))
head(boot_tci)

#####################################################################
# Precision of estimate ---------------------------------------------
#####################################################################

# Here, we are going to evaluate precise our sample estimate of the critical 
# thermal limit is, across a range of sample sizes. 
# We will use the width of the CI as a measure of precision. 

# Add width of CI to boostrapped data
ci_width_bs <- boot_tci %>%
  mutate(w_ci = upper_ci - lower_ci,
         dw_ci = abs(w_ci-lag(w_ci))) %>%
  group_by(insect_sp, sample_size) %>%
  summarise(w_ci = mean(w_ci),
           dw_ci = mean(dw_ci))
head(ci_width_bs)

# Line graph of 95% CI width
ggplot(data = ci_width_bs, aes(x = sample_size, 
                               y = w_ci,
                               colour = insect_sp)) +
  scale_colour_viridis(discrete = TRUE,
                       option = "A") +
  geom_line() + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), 
                     limits = c(0, 10)) +
  labs(x = "Sample size (n)",
       y = "Width of 95% CI (°C)") 

# Save figure to file 
ggsave("./figures/ci_ctmin_ci_precision_width.png", 
       width = 6, 
       height = 6)

#####################################################################
# Precision of estimate ---------------------------------------------
#####################################################################

# Here, we are going to evaluate precise our sample estimate of the critical 
# thermal limit is, across a range of sample sizes. 
# We will use the width of the CI as a measure of precision. 

# Line graph of 95% CI width
graph1 <- ggplot(data = ci_width_bs, aes(x = sample_size, 
                                         y = dw_ci,
                                         group = insect_sp)) +
  geom_line(colour = "grey50") + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), 
                     limits = c(0, 8)) +
  labs(x = "Sample size (n)",
       y = "Change in width of 95% CI (°C)") 
graph1

# Save figure to file 
ggsave("./figures/ci_ctmin_ci_precision_change_width.png", 
       width = 6, 
       height = 6)

# Line graph of 95% CI width (more focused y axis range)
graph2 <- ggplot(data = ci_width_bs, aes(x = sample_size, 
                                         y = dw_ci,
                                         group = insect_sp)) +
  geom_line(colour = "grey50") + 
  scale_y_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0), 
                     limits = c(0, 2)) +
  labs(x = "Sample size (n)",
       y = "Change in width of 95% CI (°C)") 
graph2

# Save figure to file 
ggsave("./figures/ci_ctmin_ci_precision_change_width_focus.png", 
       width = 6, 
       height = 6)

#####################################################################
# Precision of estimate - averaged across species -------------------
#####################################################################

# Calculate mean +- 1se change in precision 
ci_width_calc_1 <- boot_tci %>%
  group_by(insect_sp, sample_size) %>%
  mutate(w_ci = upper_ci - lower_ci) %>%
  mutate(w_ci = abs(w_ci)) %>%
  mutate(dw_ci_mean = mean(w_ci),
         dw_ci_sd = sd(w_ci),
         dw_ci_se = dw_ci_sd/sqrt(sample_size)) %>%
  slice(1) %>%
  select(insect_sp, sample_size,dw_ci_mean, dw_ci_sd, dw_ci_se) %>%
  ungroup() %>%
  group_by(sample_size) %>%
  summarise(mean_p = mean(dw_ci_mean),
            mean_se = mean(dw_ci_se))
ci_width_calc_1

# Line graph of 95% CI width
ci_width_calc_1 %>%
  dplyr::mutate(lower = mean_p - mean_se,
                upper = mean_p + mean_se) %>% 
  ggplot(data = ., aes(x = sample_size, 
                       y = mean_p)) +
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper),
              fill = "gray80") +
  geom_line() + 
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), 
                     limits = c(0, 6)) +
  labs(x = "Sample size (n)",
       y = "Width of 95% CI (°C)") 

# Save figure to file 
ggsave("./figures/ci_ctmin_ci_precision_summary.png", 
       width = 6, 
       height = 6)

#####################################################################
# Change in precision of estimate - averaged across species ---------
#####################################################################

# Calculate change in width of CI, averaged across all species 
ci_width_change <- boot_tci %>%
  dplyr::mutate(w_ci = upper_ci - lower_ci,
                ci_change = abs(w_ci-lag(w_ci))) %>%
  dplyr::group_by(insect_sp, sample_size) %>%
  dplyr::summarise(ci_change_mean = mean(ci_change, na.rm = TRUE),
                   ci_change_sd   = sd(w_ci, na.rm = TRUE),
                   ci_change_se   = ci_change_mean/sqrt(sample_size),
                   lower_se       = mean(ci_change_mean - ci_change_se),
                   upper_se       = mean(ci_change_mean + ci_change_se)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample_size) %>%
  dplyr::summarise(ci_change_mean = mean(ci_change_mean),
                   lower_se       = mean(lower_se),
                   upper_se       = mean(upper_se))
head(ci_width_change)

# Plot relationship - averaged across species
ci_width_change %>%
  ggplot(data = ., aes(x = sample_size, 
                       y = ci_change_mean)) +
  geom_ribbon(aes(ymin = lower_se, 
                  ymax = upper_se),
              fill = "gray80") +
  geom_line() + 
  #scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), 
  #                   limits = c(0, 6)) +
  labs(x = "Sample size (n)",
       y = "Change in width of 95% CI (°C)") 

# Save figure to file 
ggsave("./figures/ci_ctmin_ci_ci_width_all_species.png", 
       width = 6, 
       height = 6)

# Plot relationship - averaged across species (focused y axis)
ci_width_change %>%
  # Set upper se upper limit to 2 for ease of visualisation 
  dplyr::mutate(upper_se = dplyr::if_else(upper_se > 2, 
                                          2, 
                                          upper_se)) %>%
  dplyr::mutate(ci_change_mean = dplyr::if_else(ci_change_mean > 2, 
                                          2, 
                                          ci_change_mean)) %>%
  # Remove n = 3 for visualisation (mean > 2, so it removes the entire line)
  dplyr::filter(sample_size > 3) %>%
  ggplot(data = ., aes(x = sample_size, 
                       y = ci_change_mean)) +
  geom_ribbon(aes(ymin = lower_se, 
                  ymax = upper_se),
              fill = "gray80") +
  geom_line() + 
  scale_y_continuous(breaks = seq(0, 2, 0.25), 
                     limits = c(0, 2)) +
  labs(x = "Sample size (n)",
       y = "Change in width of 95% CI (°C)") 

# Save figure to file 
ggsave("./figures/ci_ctmin_ci_ci_width_all_species_focus.png", 
       width = 6, 
       height = 6)
