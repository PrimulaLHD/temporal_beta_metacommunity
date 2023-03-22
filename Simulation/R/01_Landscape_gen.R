
## generate the landscape location
setwd("/Users/haidong/dropbox/GitHub/temporal_beta_metacommunity")
library(mcomsimr)

landscape <- list()

set.seed(123)
for (i in 1:30) { # generating 30 landscapes, and there are 150 patches for each
  landscape[[i]] <- landscape_generate(patches = 150)
}
length(landscape)

landscape_stage_1 <- list()
for (i in 1:30) {
  landscape_stage_1[[i]] <- landscape[[i]][1:100, ]# we used the first hundred to simulate the first stage for selecting species
}

landscape_stage_2 <- list()
for (i in 1:30) {
  landscape_stage_2[[i]] <- landscape[[i]][21:120,] # we assume the first 20 patches are vanished, and new created 20 patches, namely 101-120
}

## generate environment datafram
env_df <- list()
set.seed(123)
for (i in 1:30) { 
  env_df[[i]] <- env_generate(landscape = landscape[[i]], timesteps = 2200)
}


env_df_stage_1 <- list()
for (i in 1:30) {
  env_df_stage_1[[i]] <- env_df[[i]] %>% 
    filter(patch %in% c(1:100))# we used the first hundred to simulate the first stage for selecting species
}

env_df_stage_2 <- list()
for (i in 1:30) {
  env_df_stage_2[[i]] <- env_df[[i]] %>% 
    filter(patch %in% c(21:120)) # we assume the first 20 patches are vanished, and new created 20 patches, namely 101-120
}


save(landscape,landscape_stage_1, landscape_stage_2,
     env_df, env_df_stage_1, env_df_stage_2,
     file = "Data/data2/landscape_and_env.Rdata")

load("Data/data2/landscape_and_env.Rdata")
