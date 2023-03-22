#Simulation stage 1
## by Hai-Dong Li, 2023-02-22

#############################################################################################
## part I: simulation on equal competition and narrow/flat abiotic niche


## 1. load libraries
#devtools::install_github("plthompson/mcomsimr")
library(mcomsimr)
library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

## 2. set up working directory
setwd("/Users/haidong/dropbox/GitHub/temporal_beta_metacommunity")
load("Data/data2/landscape_and_env.Rdata")

## 3. generate dispersal gradients in a log space
dispersal_gradient <- exp(seq(log(1e-5), log(1), length.out = 16))[1:15]

## 4. basic function

mcsim_eq <- function(M = 100, S = 50, dispersal, env_df, optima, rep, landscape,...){
  # env_df and landscape showed be a list object, length = number of patches
  time_sample <- seq(20, 1200, length.out = 60)
  library(dplyr)
  mcsim <- mcomsimr::simulate_MC(patches = M, species = S, dispersal = dispersal,
                                 torus = T, 
                                 env_optima = optima, 
                                 intra = 1,
                                 min_inter = 1, 
                                 max_inter = 1,
                                 env.df = env_df[[rep]],
                                 landscape = landscape[[rep]], ...)
  
  res <- mcsim$dynamics.df %>% 
    filter(time %in% time_sample) %>%
    mutate(dispersal = dispersal, reps = rep)
  res
}



## 5.simulation with paralleling
### 5.1 eq_narrow
t1 <- Sys.time()
t1
# pb <- txtProgressBar(min = 0,
#                      max = 30,
#                      style = 3,
#                      width = 15,
#                      char = "=")

for (i in 1:30) {
  cl <- makeCluster(35)
  registerDoParallel(cl)
  #(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...)
  {meta_df_eq_narrow <-foreach(dispersal = dispersal_gradient, 
                               rep = rep(i, 15),.combine='rbind') %dopar% {
                                 mcsim_eq(M = 100, S = 50, dispersal, 
                                          env_df = env_df_stage_1, rep,
                                          landscape = landscape_stage_1, 
                                          optima = narrow_optima$optima[1:50], 
                                          env_niche_breadth = 0.5)
                               };}
  
  stopCluster(cl)
  save(meta_df_eq_narrow, file = paste("Data/data2/", "meta_df_eq_narrow_", i, 
                                       ".Rdata", sep = ""))
  print(paste(i,"/30"))
  rm(meta_df_eq_narrow)
  #setTxtProgressBar(pb, rep)
}
t2 <- Sys.time()
t2
t2-t1
#close(pb)

### 5.2 eq_flat
t1 <- Sys.time()
t1
# pb <- txtProgressBar(min = 0,
#                      max = 30,
#                      style = 3,
#                      width = 15,
#                      char = "=")

for (i in 1:30) {
  cl <- makeCluster(35)
  registerDoParallel(cl)
  #(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...)
  {meta_df_eq_flat <-foreach(dispersal = dispersal_gradient, 
                               rep = rep(i, 15),.combine = 'rbind') %dopar% {
                                 mcsim_eq(M = 100, S = 50, dispersal, 
                                          env_df = env_df_stage_1, rep,
                                          landscape = landscape_stage_1, 
                                          optima = flat_optima$optima[1:50], 
                                          env_niche_breadth = 10)
                               };}
  
  stopCluster(cl)
  save(meta_df_eq_flat, file = paste("Data/data2/", "meta_df_eq_flat_", i, ".Rdata", sep = ""))
  print(paste(i,"/30"))
  rm(meta_df_eq_flat)
  #setTxtProgressBar(pb, rep)
}
t2 <- Sys.time()
t2
t2-t1
#close(pb)



#####################################################################################################################
#####################################################################################################################
# part II: simulation on stabilizing competition and narrow/flat abiotic niche
## 

## 1. load libraries
#devtools::install_github("plthompson/mcomsimr")
library(mcomsimr)
library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

## 2. set up working directory
setwd("D:/metacommunity")
load("landscape_and_env.Rdata")

## 3. generate dispersal gradients in a log space
dispersal_gradient <- exp(seq(log(1e-5), log(1), length.out = 16))[1:15]

## 4. basic function
mcsim_stab <- function(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...){
  time_sample <- seq(20, 1200, length.out = 60)
  library(dplyr)
  mcsim <- mcomsimr::simulate_MC(patches = M, species = S, dispersal = dispersal,
                                 torus = T, env.df =  env_df[[rep]], 
                                 env_optima = optima, 
                                 landscape = landscape[[rep]],...)
  
  res <- mcsim$dynamics.df %>% 
    filter(time %in% time_sample) %>%
    mutate(dispersal = dispersal, reps = rep)
  res
}

################################################################################
## 5.simulation with paralleling
################################################################################
### 5.1 stab_narrow
t1 <- Sys.time()
t1
# pb <- txtProgressBar(min = 0,
#                      max = 30,
#                      style = 3,
#                      width = 15,
#                      char = "=")

for (i in 1:30) {
  cl <- makeCluster(35)
  registerDoParallel(cl)
  #(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...)
  {meta_df_stab_narrow <-foreach(dispersal = dispersal_gradient, 
                                 rep = rep(i, 15),.combine='rbind') %dopar% {
                                   mcsim_stab(M = 100, S = 50, dispersal, 
                                              env_df = env_df_stage_1, rep,
                                              landscape = landscape_stage_1, 
                                              optima = narrow_optima$optima[1:50],
                                              int_mat = stab_int[1:50, 1:50],
                                              env_niche_breadth = 0.5)
                                 };}
  
  stopCluster(cl)
  save(meta_df_stab_narrow, file = paste("meta_df_stab_narrow_", i, 
                                         ".Rdata", sep = ""))
  print(paste(i,"/30 of stab_narrow simulation"))
  rm(meta_df_stab_narrow)
  #setTxtProgressBar(pb, rep)
}
t2 <- Sys.time()
t2
t2-t1
#close(pb)

################################################################################
### 5.2 stab_flat
t1 <- Sys.time()
t1
# pb <- txtProgressBar(min = 0,
#                      max = 30,
#                      style = 3,
#                      width = 15,
#                      char = "=")

for (i in 1:30) {
  cl <- makeCluster(35)
  registerDoParallel(cl)
  #(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...)
  {meta_df_stab_flat <-foreach(dispersal = dispersal_gradient, 
                               rep = rep(i, 15),.combine='rbind') %dopar% {
                                 mcsim_stab(M = 100, S = 50, dispersal, 
                                            env_df = env_df_stage_1, rep,
                                            landscape = landscape_stage_1, 
                                            optima = flat_optima$optima[1:50],
                                            int_mat = stab_int[1:50, 1:50],
                                            env_niche_breadth = 10)
                               };}
  
  stopCluster(cl)
  save(meta_df_stab_flat, file = paste("meta_df_stab_flat_", i, 
                                       ".Rdata", sep = ""))
  print(paste(i,"/30 of stab_flat simulation"))
  rm(meta_df_stab_flat)
  #setTxtProgressBar(pb, rep)
}
t2 <- Sys.time()
t2
t2-t1
#close(pb)






