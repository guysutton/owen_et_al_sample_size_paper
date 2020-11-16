######################################################################
######################################################################
######################################################################
# - Thermal limit studies and sample size 
# - Script 001:  Draw bootstrap resamples and clean data
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
               data.table)

# Insert instructions on how to install ThermalSampleR package
library(ThermalSampleR)

######################################################################
# Load raw data ------------------------------------------------------
######################################################################

# Import data 
sp_data <- readr::read_csv2("./data_raw/ct_min_data_raw.csv")
head(sp_data)

# Make all columns containing CT data into numeric class
sp_data_clean <- sp_data %>%
  # Clean columns names first
  janitor::clean_names() %>%
  # Make any chr columns > numeric
  dplyr::mutate(across(echthrodesis_lamorali:chrysomya_megachephala, as.numeric)) %>%
  # Remove this column
  dplyr::select(-c(chrysomya_megachephala)) %>%
  # Reshape data from wide to long format 
  tidyr::pivot_longer(
    # cols = which columns do we want to pivot/move
    cols = -rep,
    # names_to = new column name that the names of cols above will be
    # moved to. This effectively creates your categorical
    # factor levels
    names_to = "insect_sp",
    # values_to = new column where the row values of cols will be stored
    values_to = "ct_val") %>%
  # Make insect species into a factor
  dplyr::mutate(insect_sp = as.factor(insect_sp), 
                ct_val = as.numeric(ct_val)) %>%
  # Arrange new data by species
  dplyr::arrange(insect_sp) %>%
  # Filter only the first 30 reps and not NA values
  dplyr::filter(rep == c(1:30)) %>%
  tidyr::drop_na()
head(sp_data_clean)

# Because some of the data were inputed in numeric order, 
# let's randomise order of input
sp_data_clean <- sp_data_clean %>% 
  dplyr::group_by(insect_sp) %>%
  dplyr::slice(., sample(1:n())) %>%
  dplyr::mutate(sample_size = c(1:n())) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ct_val = ct_val) %>%
  dplyr::select(insect_sp, sample_size, ct_val) %>%
  dplyr::mutate(sample_size = as.numeric(sample_size), 
                insect_sp = as.factor(insect_sp))
head(sp_data_clean)

# Save raw data to disk as file called 'ct_min_data_clean.csv'
data.table::fwrite(sp_data_clean, 
                   "./data_clean/ct_min_data_clean.csv")

######################################################################
# Load cleaned data --------------------------------------------------
######################################################################

# Open the cleaned .csv file 
sp_data_clean <- readr::read_csv("./data_clean/ct_min_data_clean.csv")
head(sp_data_clean)

# Change insect species to factor
sp_data_clean <- sp_data_clean %>%
  dplyr::mutate(insect_sp = as.factor(insect_sp))
head(sp_data_clean)

#####################################################################
# Bootstrap resampling ----------------------------------------------
#####################################################################

# Will have to do bootstrap resampling in batches as my PC is 
# really slow and rubbish... 

# Extract levels of data frame to filter with 
sp_data_clean %>% 
  sapply(levels)

# First 10 groups - set.seed to make resamples reproducible 
set.seed(2012)
boot_raw_1 <- sp_data_clean %>% 
  dplyr::filter(insect_sp %in% c("Catorhintha schaffneri_APM", 
                                 "Catorhintha schaffneri_NPM", 
                                 "Dichrorampha odorata", 
                                 "Eccritotarsus catarinensis_ETG16_Summer", 
                                 "Eccritotarsus catarinensis_ETG16_Winter", 
                                 "Eccritotarsus catarinensis_ETG17", 
                                 "Eccritotarsus catarinensis_KJP_Female",  
                                 "Eccritotarsus catarinensis_KJP_Male", 
                                 "Eccritotarsus catarinensis_KTG16_Summer")) %>%
  dplyr::group_by(insect_sp) %>%
  tidyr::nest() %>%
  tidyr::crossing(sample_size = c(3:50), # No. of samples to draw per iteration
                  iter = seq(1:499)) %>% # No. of iterations 
  # Added sampling with replacement to code below
  dplyr::mutate(sample_data = purrr::map2(data, 
                                          sample_size, 
                                          ~ dplyr::sample_n(.x, .y, 
                                                            replace = TRUE))) %>% 
  dplyr::mutate(calc = purrr::map(sample_data, 
                                  ~dplyr::summarise(., 
                                                    mean = mean(ct_val), 
                                                    sd = sd((ct_val))))) %>% 
  dplyr::select(insect_sp, sample_size, iter, calc) %>% 
  tidyr::unnest(cols = calc)

# Check the bootstrapped dataset
head(boot_raw_1)

# Second 10 groups 
set.seed(2012)
boot_raw_2 <- sp_data_clean %>% 
  dplyr::filter(insect_sp %in% c("Eccritotarsus catarinensis_KTG16_Winter",
                                 "Eccritotarsus catarinensis_KTG17", 
                                 "Eccritotarsus catarinensis_QJP_Female", 
                                 "Eccritotarsus catarinensis_QTG", 
                                 "Eccritotarsus catarinensis_WJC", 
                                 "Eccritotarsus catarinensisQJP Male", 
                                 "Eccritotarsus eichhorniae_QTG", 
                                 "Eccritotarsus eichhorniae_WCO", 
                                 "Echthrodesis lamorali", 
                                 "Hydrellia egeriae_RS_Female")) %>%
  dplyr::group_by(insect_sp) %>%
  tidyr::nest() %>%
  tidyr::crossing(sample_size = c(3:50), # No. of samples to draw per iteration
                  iter = seq(1:499)) %>% # No. of iterations 
  # Added sampling with replacement to code below
  dplyr::mutate(sample_data = map2(data, 
                                   sample_size, 
                                   ~ dplyr::sample_n(.x, .y, 
                                                     replace = TRUE))) %>% 
  dplyr::mutate(calc = purrr::map(sample_data, 
                                  ~dplyr::summarise(., 
                                                    mean = mean(ct_val), 
                                                    sd = sd((ct_val))))) %>% 
  dplyr::select(insect_sp, sample_size, iter, calc) %>% 
  tidyr::unnest(cols = calc)

# Check the bootstrapped dataset
head(boot_raw_2)

# Third 10 groups 
set.seed(2012)
boot_raw_3 <- sp_data_clean %>% 
  dplyr::filter(insect_sp %in% c("Hydrellia egeriae_RS_Male",
                                 "Megamelus scutellaris_DCT_Female",
                                 "Megamelus scutellaris_DCT_Male",
                                 "Megamelus scutellaris_KCT_Female",
                                 "Megamelus scutellaris_KCT_Male",
                                 "Megamelus scutellaris_WCT_Female", 
                                 "Megamelus scutellaris_WCT_Male",
                                 "Megamelus scuttelaris_QJC",
                                 "Neochetina bruchi",
                                 "Neochetina eichhorniae_ECO_Female")) %>%
  dplyr::group_by(insect_sp) %>%
  tidyr::nest() %>%
  tidyr::crossing(sample_size = c(3:50), # No. of samples to draw per iteration
                  iter = seq(1:499)) %>% # No. of iterations 
  # Added sampling with replacement to code below
  dplyr::mutate(sample_data = purrr::map2(data, 
                                          sample_size, 
                                          ~ dplyr::sample_n(.x, .y, 
                                                            replace = TRUE))) %>% 
  dplyr::mutate(calc = purrr::map(sample_data, 
                                  ~ dplyr::summarise(., 
                                                     mean = mean(ct_val),
                                                     sd = sd((ct_val))))) %>% 
  dplyr::select(insect_sp, sample_size, iter, calc) %>% 
  tidyr::unnest(cols = calc)

# Check the bootstrapped dataset
head(boot_raw_3)

# Fourth group 
set.seed(2012)
boot_raw_4 <- sp_data_clean %>% 
  dplyr::filter(insect_sp %in% c("Neochetina eichhorniae_ECO_Male", 
                                 "Neochetina eichhorniae_KCO_Female", 
                                 "Neochetina eichhorniae_KCO_Male", 
                                 "Neochetina eichhorniae_KDR_Female",   
                                 "Neochetina eichhorniae_SDR_Female", 
                                 "Neochetina eichhorniae_SDR_Male", 
                                 "Neochetina eichhorniae_WJC", 
                                 "Neochetina eichhorniaeKDR Males",  
                                 "Neohydronomus affinis", 
                                 "Stenopelmus rufinasus",             
                                 "Thaumatotibia leucotreta_std_diet_CO_30", 
                                 "Thaumatotibia leucotreta_tre_diet_CO_30")) %>%
  dplyr::group_by(insect_sp) %>%
  tidyr::nest() %>%
  tidyr::crossing(sample_size = c(3:50), # No. of samples to draw per iteration 
                  iter = seq(1:499)) %>% # No. of iterations 
  # Added sampling with replacement to code below
  dplyr::mutate(sample_data = purrr::map2(data, 
                                          sample_size, 
                                          ~ dplyr::sample_n(.x, .y, 
                                                            replace = TRUE))) %>% 
  dplyr::mutate(calc = purrr::map(sample_data, 
                                  ~ dplyr::summarise(., 
                                                     mean = mean(ct_val), 
                                                     sd = sd((ct_val))))) %>% 
  dplyr::select(insect_sp, sample_size, iter, calc) %>% 
  tidyr::unnest(cols = calc)

# Check the bootstrapped dataset
head(boot_raw_4)

# Combine all the bootstrapped datasets into one file 
boot_raw <- dplyr::bind_rows(boot_raw_1, 
                             boot_raw_2, 
                             boot_raw_3, 
                             boot_raw_4)
head(boot_raw)

# Save raw bootstrap samples to disk 
data.table::fwrite(boot_raw, 
                   "./data_clean/ct_min_bootstrap_samples_raw.csv")

#####################################################################
# Load bootstrap samples --------------------------------------------
#####################################################################

# Load raw bootstrap data 
boot_raw <- readr::read_csv("./data_clean/ct_min_bootstrap_samples_raw.csv")
head(boot_raw)

# Change insect species to factor
boot_raw <- boot_raw %>%
  dplyr::mutate(insect_sp = as.factor(insect_sp))
head(boot_raw)

# Remove Hydrellia egeriae as they were preliminary data 
# and Chrysomya megacephala (not sure about that data)
boot_raw <- boot_raw %>%
  dplyr::filter(!insect_sp %in% c("Hydrellia egeriae_RS_Female",
                                  "Hydrellia egeriae_RS_Male",
                                  "Chrysomya megachephala"))

######################################################################
# Find population parameters -----------------------------------------
######################################################################

# Get median values per species (median for n = 30; 
# max. sample size in our study)
# Following Pearson and Groves (2013), we assumed that our maximum 
# sample size (n = 30) is a reasonable approximation of the actual 
# critical thermal limit value (population parameter)
# i.e. if we hypothetically had sampled the entire population
# Obviously, we must interpret these data with caution, as there is no 
# possible way to 100% accurately determine the population parameter, 
# so we must derive an ESTIMATE.

# Estimate population CT value per species 
median_vals <- boot_raw %>%
  dplyr::filter(sample_size == 30) %>%
  dplyr::group_by(insect_sp) %>%
  dplyr::mutate(median_pop_val = median(mean)) %>%
  dplyr::slice(1) %>%
  dplyr::select(insect_sp, median_pop_val)

# Add median value per species to the bootstrapped dataset
# Join the two datsets by the 'insect_sp' column
boot_proc <- dplyr::full_join(boot_raw, 
                              median_vals, 
                              by = "insect_sp") 

######################################################################
# Calculate confidence intervals  ------------------------------------
######################################################################

# Because our sample sizes are < 30, we cannot assume the central limit 
# theorem to calculate standard 95% confidence intervals (hereafter 'CI') 
# (e.g. from a normal distribution), so we will use the t-distribution 
# to calculate 95% CI's. 

# We use Student's t-distribution to calculate a 95% CI: 
# General formula: mean +- t(alpha)/2 x s.d./sqrt(n);
# So, we need to specify a different t-value for each sample size, 
# to accurately reflect the degrees of freedom (here = n - 1)

boot_tci <- boot_proc %>%
  # Add standard errors
  dplyr::mutate(std_error = sd/sqrt(sample_size)) %>%
  # Now find error 
  dplyr::mutate(error = qt(0.975, df = sample_size - 1) * std_error) %>%
  # Calculate lower and upper 95% CI limits
  dplyr::mutate(lower_ci = mean - error,
                upper_ci = mean + error) %>%
  dplyr::ungroup()

# Check code worked 
head(boot_tci)

# Save this file to disc
data.table::fwrite(boot_tci, 
                   "./data_clean/ct_min_bootstrap_with_ci.csv")
