library(tidyverse)
library(RColorBrewer)
library(glmmTMB)

# The occurence data is stored somewhere and loaded:
# alldata <- read.table("occurrence.txt", sep="\t",header=T, fileEncoding = "UTF8")

# Data cleaning:
source("data_cleaning.R")

data_agg <- read_rds("data_agg_depth.rds")

data_agg <- data_agg %>% 
  mutate(year = year - min(year))

options(max.print = 1e4)

data_agg %>% 
  ggplot(aes(x = time, y = log(individualCount))) +
  # geom_line(aes(group = interaction(species, year, locality, depth), colour = depth)) +
  # geom_point(aes(colour = depth)) +
  geom_line(aes(group = interaction(species, year, locality, depth), colour = locality)) +
  geom_point(aes(colour = locality)) +
  facet_grid(phylum + species ~ year)

# If the figure is too big to make sense of it, you can plot smaller sections of it, just
# un-comment one of the filters:
data_agg %>% 
  # filter(year < 13,
  #        phylum == "Rotifera") %>% 
  # filter(year > 12,
  #        phylum == "Rotifera") %>% 
  # filter(year < 13,
  #        phylum != "Rotifera") %>%
  # filter(year > 12,
  #        phylum != "Rotifera") %>%
  ggplot(aes(x = time, y = log(individualCount))) +
  # geom_line(aes(group = interaction(species, year, locality, depth), colour = depth)) +
  # geom_point(aes(colour = depth)) +
  geom_line(aes(group = interaction(species, year, locality, depth), colour = locality)) +
  geom_point(aes(colour = locality)) +
  facet_grid(phylum + species ~ year)

options(max.print = 1e4)

# Depending on the data set, we specify different models that take into account
# the number localities and phylum:
f_est_m1 <- function(data, ver = "b"){
  # ver = "a" is for data set that contains more than one locality and phylum:
  if (ver == "a"){
    mod <- glmmTMB(individualCount ~ 
                     # species-specific response to annual environmental
                     # fluctuations within locations:
                     ou(nf_year + 0 | species:locality) +
                     # general response to annual environmental fluctuations:
                     ou(nf_year + 0 | lake) +
                     # species heterogeneity:
                     (1 | species) +
                     # local species heterogeneity:
                     (1 | species:locality) +
                     # fixed effect difference between phylum
                     phylum +
                     # observation-level random effect:
                     (1 | species:locality:time:nf_year:depth),
                   # observation model
                   family = truncated_poisson(link = "log"),
                   data = add_column(data, lake = "A"))
  }
  # ver = "b" is for data set that contains more than one locality, and a single phylum
  else if (ver == "b"){
    mod <- glmmTMB(individualCount ~ 
                     # species-specific response to annual environmental
                     # fluctuations within locations:
                     ou(nf_year + 0 | species:locality) +
                     # general response to annual environmental fluctuations:
                     ou(nf_year + 0 | lake) +
                     # species heterogeneity:
                     (1 | species) +
                     # local species heterogeneity
                     (1 | species:locality) +
                     # observation-level random effect:
                     (1 | species:locality:time:nf_year:depth),
                   family = truncated_poisson(link = "log"),
                   data = add_column(data, lake = "A"))
  }
  # ver = "d" is for data set that contains one location and two phylum
  else if (ver == "d"){
    mod <- glmmTMB(individualCount ~ 
                     # species-specific response to annual environmental
                     # fluctuations within locations:
                     ou(nf_year + 0 | species:locality) +
                     # general response to annual environmental fluctuations
                     ou(nf_year + 0 | lake) +
                     # species heterogeneity:
                     (1 | species) +
                     # fixed effect difference between phylum
                     phylum +
                     # observation-level random effect:
                     (1 | species:locality:time:nf_year:depth),
                   # observation model:
                   family = truncated_poisson(link = "log"),
                   data = add_column(data, lake = "A"))
  }
  # ver = "c" is for data set at one location and a single phylum:
  else{
    mod <- glmmTMB(individualCount ~ 
                     # species-specific response to annual environmental
                     # fluctuations within locations:
                     ou(nf_year + 0 | species:locality) +
                     # general response to annual environmental fluctuations:
                     ou(nf_year + 0 | lake) +
                     # species heterogeneity:
                     (1 | species) +
                     # observation-level random effect:
                     (1 | species:locality:time:nf_year:depth),
                   # observation model:
                   family = truncated_poisson(link = "log"),
                   data = add_column(data, lake = "A"))
  }
  
  return(mod)
}

# Again, depending on the data set, we make a function that extract the different
# parameters:
# regardless of which "ver" is selected, you need to specify the model and in
# case of ver = "b" and "c" you need to specify which phylum is in the model.
# ver = "a" is for data set that contains more than one locality and phylum:
# ver = "b" is for data set that contains more than one locality, and a single phylum
# ver = "d" is for data set that contains one location and two phylum
# ver = "c" is for data set at one location and a single phylum:
f_par_m1 <- function(mod, ver = "b", phyl = "arthropoda"){
  if (ver == "a"){
    tmp <- mod$fit$par[1:2]
    tmp <- c(tmp, exp(mod$fit$par[-c(1:2)])^c(rep(c(2, -1), 2), rep(2, 3)))
    names(tmp) <- c("beta_arthropoda", # Fixed effect for arthropoda phylum (mean log abundance)
                    "beta_rotifera", # Difference between arthropoda and rotifera phylum
                    "year_env_sp_loc", # species-specific response to yearly 
                    # environmental fluctuations (within locations)
                    "tr_ysl", # mean return to equilibrium (the inverse of the 
                    # strength of density regulation) for the species-specific 
                    # response to yearly environmental fluctuations (within locations)
                    "year_env_lake", # general response to yearly 
                    # environmental fluctuations across all locations
                    "tr_lake", # mean return to equilibrium (the inverse of the 
                    # strength of density regulation) for the general response 
                    # to yearly environmental fluctuations across all locations
                    "het", # species heterogeneity (across all location)
                    "het_loc", # local species heterogeneity
                    "obs") # observation-level random effect (variation in 
                    # observations within a year)
  }
  else if (ver == "b"){
    tmp <- mod$fit$par[1]
    tmp <- c(tmp, exp(mod$fit$par[-c(1)])^c(rep(c(2, -1), 2), rep(2, 3)))
    names(tmp) <- c(paste("beta", phyl, sep = "_"), # fixed effect for specified phylum
                    "year_env_sp_loc", "tr_ysl", # see above
                    "year_env_lake", "tr_lake",
                    "het", "het_loc", "obs")
  }
  else if (ver == "d"){
    tmp <- mod$fit$par[1:2]
    tmp <- c(tmp, exp(mod$fit$par[-c(1:2)])^c(rep(c(2, -1), 2), rep(2, 2)))
    names(tmp) <- c("beta_arthropoda", "beta_rotifera", # see above
                    "year_env_sp_loc", "tr_ysl",
                    "year_env_lake", "tr_lake",
                    "het", "obs")
  }
  else {
    tmp <- mod$fit$par[1]
    tmp <- c(tmp, exp(mod$fit$par[-c(1)])^c(rep(c(2, -1), 2), rep(2, 2)))
             names(tmp) <- c(paste("beta", phyl, sep = "_"), # see above
                             "year_env_sp_loc", "tr_ysl",
                             "year_env_lake", "tr_lake",
                             "het", "obs")
  }
  return(tmp)
}

# Fitting all the different models, although we will only use "a" and "d" for now.
# All data:
mod1a <- f_est_m1(data_agg, ver = "a")
# Single phylum:
mod1b1 <- f_est_m1(filter(data_agg, phylum == unique(data_agg$phylum)[1]), ver = "b")
mod1b2 <- f_est_m1(filter(data_agg, phylum == unique(data_agg$phylum)[2]), ver = "b")
# Single phylum, single locality:
mod1c11 <- f_est_m1(filter(data_agg, 
                           phylum == unique(data_agg$phylum)[1], 
                           locality == unique(data_agg$locality)[1]), ver = "c")
mod1c12 <- f_est_m1(filter(data_agg, 
                           phylum == unique(data_agg$phylum)[1], 
                           locality == unique(data_agg$locality)[2]), ver = "c")
mod1c13 <- f_est_m1(filter(data_agg, 
                           phylum == unique(data_agg$phylum)[1], 
                           locality == unique(data_agg$locality)[3]), ver = "c")
mod1c21 <- f_est_m1(filter(data_agg, 
                           phylum == unique(data_agg$phylum)[2], 
                           locality == unique(data_agg$locality)[1]), ver = "c")
mod1c22 <- f_est_m1(filter(data_agg, 
                           phylum == unique(data_agg$phylum)[2], 
                           locality == unique(data_agg$locality)[2]), ver = "c")
mod1c23 <- f_est_m1(filter(data_agg, 
                           phylum == unique(data_agg$phylum)[2], 
                           locality == unique(data_agg$locality)[3]), ver = "c")
# Single locality:
mod1d1 <- f_est_m1(filter(data_agg, 
                          locality == unique(data_agg$locality)[1]), ver = "d")
mod1d2 <- f_est_m1(filter(data_agg, 
                          locality == unique(data_agg$locality)[2]), ver = "d")
mod1d3 <- f_est_m1(filter(data_agg, 
                          locality == unique(data_agg$locality)[3]), ver = "d")

# Extract parameters and put everything in a table:
tmp_tbl <- bind_rows(
  add_column(as_tibble_row(f_par_m1(mod1a, ver = "a")), phylum = "Both", locality = "All"),
  add_column(as_tibble_row(f_par_m1(mod1b1, ver = "b")), phylum = "Arthopoda", locality = "All"),
  add_column(as_tibble_row(f_par_m1(mod1b2, ver = "b", phyl = "rotifera")), phylum = "Rotifera", locality = "All"),
  add_column(as_tibble_row(f_par_m1(mod1c11, ver = "c")), phylum = "Arthopoda", locality = "Kilvatn"),
  add_column(as_tibble_row(f_par_m1(mod1c13, ver = "c")), phylum = "Arthopoda", locality = "Lille Jonsvatn"),
  add_column(as_tibble_row(f_par_m1(mod1c12, ver = "c")), phylum = "Arthopoda", locality = "Store Jonsvatn"),
  add_column(as_tibble_row(f_par_m1(mod1c21, ver = "c", phyl = "rotifera")), phylum = "Rotifera", locality = "Kilvatn"),
  add_column(as_tibble_row(f_par_m1(mod1c22, ver = "c", phyl = "rotifera")), phylum = "Rotifera", locality = "Lille Jonsvatn"),
  add_column(as_tibble_row(f_par_m1(mod1c23, ver = "c", phyl = "rotifera")), phylum = "Rotifera", locality = "Store Jonsvatn"),
  add_column(as_tibble_row(f_par_m1(mod1d1, ver = "d", phyl = "Both")), phylum = "Both", locality = "Kilvatn"),
  add_column(as_tibble_row(f_par_m1(mod1d2, ver = "d", phyl = "Both")), phylum = "Both", locality = "Lille Jonsvatn"),
  add_column(as_tibble_row(f_par_m1(mod1d3, ver = "d", phyl = "Both")), phylum = "Both", locality = "Store Jonsvatn"))

tmp_tbl 

# Make a function that simulates new sample and estimate the parameter (parametric bootstrap)
f_sim_m1 <- function(data, mod, ver = "b"){
  if (ver == "a"){
    sim_est <- glmmTMB(sim_1 ~ 
                     # species-specific response to annual environmental
                     # fluctuations within locations:
                     ou(nf_year + 0 | species:locality) +
                     # general response to annual environmental fluctuations
                     ou(nf_year + 0 | lake) +
                     # species heterogeneity:
                     (1 | species) +
                     # local species heterogeneity:
                     (1 | species:locality) +
                     # fixed effect among phylum:
                     phylum +
                     # observation-level random effect:
                     (1 | species:locality:time:nf_year:depth),
                     # observation model:
                   family = truncated_poisson(link = "log"),
                   data = bind_cols(add_column(data, lake = "A"),
                                    simulate(mod)))
  }
  else if (ver == "b"){
    sim_est <- glmmTMB(sim_1 ~ 
                     # species-specific response to annual environmental
                     # fluctuations within locations:
                     ou(nf_year + 0 | species:locality) +
                     # general response to annual environmental fluctuations
                     ou(nf_year + 0 | lake) +
                     # species heterogeneity:
                     (1 | species) +
                     # local species heterogeneity:
                     (1 | species:locality) +
                     # observation-level random effect:
                     (1 | species:locality:time:nf_year:depth),
                     # observation model:
                   family = truncated_poisson(link = "log"),
                   data = bind_cols(add_column(data, lake = "A"),
                                    simulate(mod)))
  }
  else if (ver == "d"){
    sim_est <- glmmTMB(sim_1 ~ 
                         # species-specific response to annual environmental
                         # fluctuations within locations:
                         ou(nf_year + 0 | species:locality) +
                         # general response to annual environmental fluctuations
                         ou(nf_year + 0 | lake) +
                         # species heterogeneity:
                         (1 | species) +
                         # fixed effect among phylum
                         phylum +
                         # observation-level random effect:
                         (1 | species:locality:time:nf_year:depth),
                       # observation model:
                       family = truncated_poisson(link = "log"),
                       data = bind_cols(add_column(data, lake = "A"),
                                        simulate(mod)))
  }
  else{
    sim_est <- glmmTMB(sim_1 ~ 
                     # species-specific response to annual environmental
                     # fluctuations within locations:
                     ou(nf_year + 0 | species:locality) +
                     # general response to annual environmental fluctuations
                     ou(nf_year + 0 | lake) +
                     # species heterogeneity:
                     (1 | species) +
                     # observation-level random effect:
                     (1 | species:locality:time:nf_year:depth),
                     # observation model:
                   family = truncated_poisson(link = "log"),
                   data = bind_cols(add_column(data, lake = "A"),
                                    simulate(mod)))
  }
  
  return(sim_est)
}

# Test of simulation and estimation and compare to intial result:
# sim_m1a_1 <- f_sim_m1(data_agg, mod1a, ver = "a")
# f_par_m1(sim_m1a_1, ver = "a")
# f_par_m1(mod1a, ver = "a")
# sim_m1b1_1 <- f_sim_m1(filter(data_agg, phylum == unique(data_agg$phylum)[1]), mod1b1, ver = "b")
# f_par_m1(sim_m1b1_1, ver = "b")
# f_par_m1(mod1b1, ver = "b")
# sim_m1b2_1 <- f_sim_m1(filter(data_agg, phylum == unique(data_agg$phylum)[2]), mod1b2, ver = "b")
# f_par_m1(sim_m1b2_1, ver = "b")
# f_par_m1(mod1b2, ver = "b")
# sim_m1c11_1 <- f_sim_m1(filter(data_agg, 
#                                phylum == unique(data_agg$phylum)[1], 
#                                locality == unique(data_agg$locality)[1]), mod1c11, ver = "c")
# f_par_m1(sim_m1c11_1, ver = "c")
# f_par_m1(mod1c11, ver = "c")
# sim_m1c12_1 <- f_sim_m1(filter(data_agg, 
#                                phylum == unique(data_agg$phylum)[1], 
#                                locality == unique(data_agg$locality)[2]), mod1c12, ver = "c")
# f_par_m1(sim_m1c12_1, ver = "c")
# f_par_m1(mod1c12, ver = "c")
# sim_m1c13_1 <- f_sim_m1(filter(data_agg, 
#                                phylum == unique(data_agg$phylum)[1], 
#                                locality == unique(data_agg$locality)[3]), mod1c13, ver = "c")
# f_par_m1(sim_m1c13_1, ver = "c")
# f_par_m1(mod1c13, ver = "c")
# sim_m1c21_1 <- f_sim_m1(filter(data_agg, 
#                                phylum == unique(data_agg$phylum)[2], 
#                                locality == unique(data_agg$locality)[1]), mod1c21, ver = "c")
# f_par_m1(sim_m1c21_1, ver = "c")
# f_par_m1(mod1c21, ver = "c")
# sim_m1c22_1 <- f_sim_m1(filter(data_agg, 
#                                phylum == unique(data_agg$phylum)[2], 
#                                locality == unique(data_agg$locality)[2]), mod1c22, ver = "c")
# f_par_m1(sim_m1c22_1, ver = "c")
# f_par_m1(mod1c22, ver = "c")
# sim_m1c23_1 <- f_sim_m1(filter(data_agg, 
#                                phylum == unique(data_agg$phylum)[2], 
#                                locality == unique(data_agg$locality)[3]), mod1c23, ver = "c")
# f_par_m1(sim_m1c23_1, ver = "c")
# f_par_m1(mod1c23, ver = "c")

# Run things in parallel for parametric bootstrap:
library(foreach)
library(parallel)
library(doParallel)

index <- 1:1000

# Fit all models in one go:
# You can comment out models that you do not need.
f_sim_run <- function(i){
  
  sim_m1a_1 <- f_sim_m1(data_agg, mod1a, ver = "a")
  sim_m1b1_1 <- f_sim_m1(filter(data_agg, phylum == unique(data_agg$phylum)[1]), mod1b1, ver = "b")
  sim_m1b2_1 <- f_sim_m1(filter(data_agg, phylum == unique(data_agg$phylum)[2]), mod1b2, ver = "b")
  sim_m1c11_1 <- f_sim_m1(filter(data_agg, 
                                 phylum == unique(data_agg$phylum)[1], 
                                 locality == unique(data_agg$locality)[1]), mod1c11, ver = "c")
  sim_m1c12_1 <- f_sim_m1(filter(data_agg, 
                                 phylum == unique(data_agg$phylum)[1], 
                                 locality == unique(data_agg$locality)[2]), mod1c12, ver = "c")
  sim_m1c13_1 <- f_sim_m1(filter(data_agg, 
                                 phylum == unique(data_agg$phylum)[1], 
                                 locality == unique(data_agg$locality)[3]), mod1c13, ver = "c")
  sim_m1c21_1 <- f_sim_m1(filter(data_agg, 
                                 phylum == unique(data_agg$phylum)[2], 
                                 locality == unique(data_agg$locality)[1]), mod1c21, ver = "c")
  sim_m1c22_1 <- f_sim_m1(filter(data_agg, 
                                 phylum == unique(data_agg$phylum)[2], 
                                 locality == unique(data_agg$locality)[2]), mod1c22, ver = "c")
  sim_m1c23_1 <- f_sim_m1(filter(data_agg, 
                                 phylum == unique(data_agg$phylum)[2], 
                                 locality == unique(data_agg$locality)[3]), mod1c23, ver = "c")
  sim_m1d1_1 <- f_sim_m1(filter(data_agg, 
                                locality == unique(data_agg$locality)[1]), mod1d1, ver = "d")
  sim_m1d2_1 <- f_sim_m1(filter(data_agg, 
                                locality == unique(data_agg$locality)[2]), mod1d2, ver = "d")
  sim_m1d3_1 <- f_sim_m1(filter(data_agg, 
                                locality == unique(data_agg$locality)[3]), mod1d3, ver = "d")
  
  tmp_tbl <- bind_rows(
    add_column(as_tibble_row(f_par_m1(sim_m1a_1, ver = "a")), phylum = "Both", locality = "All"),
    add_column(as_tibble_row(f_par_m1(sim_m1b1_1, ver = "b")), phylum = "Arthropoda", locality = "All"),
    add_column(as_tibble_row(f_par_m1(sim_m1b2_1, ver = "b", phyl = "rotifera")), phylum = "Rotifera", locality = "All"),
    add_column(as_tibble_row(f_par_m1(sim_m1c11_1, ver = "c")), phylum = "Arthropoda", locality = "Kilvatn"),
    add_column(as_tibble_row(f_par_m1(sim_m1c12_1, ver = "c")), phylum = "Arthropoda", locality = "Lille Jonsvatn"),
    add_column(as_tibble_row(f_par_m1(sim_m1c13_1, ver = "c")), phylum = "Arthropoda", locality = "Store Jonsvatn"),
    add_column(as_tibble_row(f_par_m1(sim_m1c21_1, ver = "c", phyl = "rotifera")), phylum = "Rotifera", locality = "Kilvatn"),
    add_column(as_tibble_row(f_par_m1(sim_m1c22_1, ver = "c", phyl = "rotifera")), phylum = "Rotifera", locality = "Lille Jonsvatn"),
    add_column(as_tibble_row(f_par_m1(sim_m1c23_1, ver = "c", phyl = "rotifera")), phylum = "Rotifera", locality = "Store Jonsvatn"),
    add_column(as_tibble_row(f_par_m1(sim_m1d1_1, ver = "d", phyl = "Both")), phylum = "Both", locality = "Kilvatn"),
    add_column(as_tibble_row(f_par_m1(sim_m1d2_1, ver = "d", phyl = "Both")), phylum = "Both", locality = "Lille Jonsvatn"),
    add_column(as_tibble_row(f_par_m1(sim_m1d3_1, ver = "d", phyl = "Both")), phylum = "Both", locality = "Store Jonsvatn"))
  
  tmp_conv <- c(sim_m1a_1$fit$convergence,
                sim_m1b1_1$fit$convergence,
                sim_m1b2_1$fit$convergence,
                sim_m1c11_1$fit$convergence,
                sim_m1c12_1$fit$convergence,
                sim_m1c13_1$fit$convergence,
                sim_m1c21_1$fit$convergence,
                sim_m1c22_1$fit$convergence,
                sim_m1c23_1$fit$convergence,
                sim_m1d1_1$fit$convergence,
                sim_m1d2_1$fit$convergence,
                sim_m1d3_1$fit$convergence)
  
  tmp_aic <- c(AIC(sim_m1a_1),
               AIC(sim_m1b1_1),
               AIC(sim_m1b2_1),
               AIC(sim_m1c11_1),
               AIC(sim_m1c12_1),
               AIC(sim_m1c13_1),
               AIC(sim_m1c21_1),
               AIC(sim_m1c22_1),
               AIC(sim_m1c23_1),
               AIC(sim_m1d1_1),
               AIC(sim_m1d2_1),
               AIC(sim_m1d3_1))
  
  tmp_tbl <- tmp_tbl %>% 
    add_column(conv = tmp_conv,
               aic = tmp_aic)
  
  return(tmp_tbl)
  
}

f_sim_run(index[1])

cores <- detectCores()
clusters <- makeCluster(cores - 1)
registerDoParallel(clusters)

run01 <- foreach(i = index,
                 .combine = "rbind",
                 .packages = c("dplyr", "tidyverse", "glmmTMB"),
                 .errorhandling = "remove") %dopar% {f_sim_run(i)}

stopCluster(clusters)

# write_rds(run01, "parallel_run_with_depth.rds")
# run01 <- read_rds("parallel_run04.rds")
run01 <- read_rds("parallel_run_with_depth.rds")
# Add bootstrap replicate index
run01 <- run01 %>%
  add_column(boot = rep(1:1000, each = 12))

# Just a glimpse of the different results:
run01 %>% 
  pivot_longer(beta_arthropoda:obs) %>% 
  filter(!is.na(value)) %>% 
  ggplot(aes(x = phylum, y = value)) +
  geom_point(aes(colour = phylum)) + 
  facet_grid(name ~ locality, scales = "free_y")

 # Calculate different summary statistics:
library(HDInterval)
library(modeest)

f_sum <- function(x, p = 0.95){
  x <- x$value
  tmp <- c(mean(x), sd(x), median(x), meanshift(x), hdi(x, p),
           quantile(x, probs = c((1 - p)/2, 1 - (1 - p) / 2)))
  names(tmp)[1:4] <- c("mean", "sd", "median", "mode")
  return(as_tibble_row(tmp))
}

#run01_sum <- run01 %>% 
#  pivot_longer(beta_arthropoda:obs) %>% 
#  filter(!is.na(value)) %>% 
#  group_by(locality, phylum, name) %>% 
#  nest() %>% 
#  mutate(stats = map(data, f_sum)) %>% 
#  unnest(stats) %>% 
#  select(-data) %>% 
#  ungroup()

library(HDInterval)
library(modeest)

run01_sum <- run01 %>% 
  filter(!(phylum == "Both" &
             locality == "Kilvatn" &
             boot %in% c("200"))) %>%
  pivot_longer(beta_arthropoda:obs) %>% 
  filter(!is.na(value)) %>% 
  group_by(locality, phylum, name) %>% 
  nest() %>% 
  mutate(stats = map(data, f_sum)) %>% 
  unnest(stats) %>% 
  select(-data) %>% 
  ungroup()


name_tbl <- bind_rows(tibble(name = c("tr_ysl", "tr_lake"),
       name_wide = c("Mean return time for species-specific response",
                     "Mean return time for general response")),
       tibble(name = c("obs", "year_env_sp_loc", "year_env_lake", "het", "het_loc"),
       name_wide = c("Observational/within-year variance",
                     "Species-specific environmental variance",
                     "Common environmental variance", 
                     "Species heterogeneity", 
                     "Local species heterogeneity")),
       tibble(name = c("beta_arthropoda", "beta_rotifera"),
                 name_wide = c("Mean log abundance for arthropoda",
                               "Difference in mean log abundance between arthropoda and rotifera")))

library(openxlsx)
xlxs_tbl <- left_join(filter(run01_sum, phylum == "Both"), name_tbl, by = join_by(name)) %>% 
  select(locality, name_wide, mean, sd, median, mode, lower, upper, `2.5%`, `97.5%`)
wb <- createWorkbook()
addWorksheet(wb, sheetName = "results_with_depth")
writeData(wb, sheet = "results_with_depth", x = xlxs_tbl)
setColWidths(wb, sheet = "results_with_depth", cols = 1:ncol(xlxs_tbl), widths = "auto")
saveWorkbook(wb, "results_with_depth.xlsx")

# Add some classes of different parameters
para_tbl <- tibble(name = unique(run01_sum$name),
                   parameter = c(rep("Fixed effect", 2), rep(c("Variance", "Dynamics"), 2), rep ("Variance", 3)))



# Finally, we plot the variance components for each of the four data sets:
# Paper fig. 5

tmp_var <- tmp_tbl %>% 
  pivot_longer(beta_arthropoda:obs) %>% 
  filter(name %in% c("year_env_sp_loc", "year_env_lake", "het", "het_loc", "obs")) %>% 
  mutate(name = factor(name, levels = c("year_env_sp_loc", "year_env_lake", "het", "het_loc", "obs")))

var_comp <- tmp_var %>%
  left_join(tibble(name = c("obs", "year_env_sp_loc", "year_env_lake", "het", "het_loc"),
                   name_wide = c("Observational/within-year variance",
                                 "Species-specific environmental variance",
                                 "Common environmental variance", 
                                 "Species heterogeneity", 
                                 "Local species heterogeneity")),
            by = join_by(name)) %>% 
  mutate(name_wide = factor(name_wide,
                            levels = c("Observational/within-year variance", 
                                       "Species-specific environmental variance", 
                                       "Common environmental variance", 
                                       "Species heterogeneity",
                                       "Local species heterogeneity"))) %>% 
  filter(phylum == "Both") %>% 
  ggplot(aes(x = locality, y = value)) + 
  geom_col(aes(fill = name_wide), colour = "black") +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(expand = expansion(add = c(0, 0.1)),
                     breaks = seq(0, 7, by = 1)) +
  guides(fill = guide_legend("Variance component:")) +
  # scale_fill_brewer(type = "qual") +
  labs(x = "",
       y = "Total variance") +
  theme_bw() +
  theme(axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5)), 
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5)),
        legend.position = "bottom",
        legend.direction = "vertical")# +
  # facet_wrap(~ locality, nrow = 1, scales = "free_x")

var_comp

#ggsave("fig_var_components_stacked.jpg", dpi = 600, width = 8000, height = 6000, units = "px")