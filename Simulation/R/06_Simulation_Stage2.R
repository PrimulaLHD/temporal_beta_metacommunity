# simulation stage 2
## function to simulate species and patch change
MC_change <- function(MC_out, dispersal_gradient, optima, int_mat){
  
  MC_t1 <- MC_out %>% 
    filter(time == 1200) %>% # select the last step
    filter(patch %in% 21:100) %>% ## remove the first 20 patches
    filter(N > 0)  # keep species that survive
  #spe_pool <- list()
  S_size <- vector(length = 15)
  spe_inf <- list()
  optima_new <- list()
  int_mat_new <- list()
  
  for (i in 1:15){
    species_left_0 <- unique((MC_t1 %>% 
                                filter(dispersal == dispersal_gradient[i]))$`species`)
    if(length(species_left_0) > 2) {
      species_left <- species_left_0[3:length(species_left_0)] # remove 2 species
    new_pool <- c(species_left, 51:70) # introduce 20 species to meta-community
    } else{
      new_pool <- c(51:70)
    }
    
    S_size[i] <- length(new_pool)
    spe_inf[[i]] <- data.frame(dispersal = dispersal_gradient[i], 
                               species = 1:length(new_pool), species2 = new_pool)
    optima_new_df <- optima %>% 
      filter(species %in% new_pool)
    
    optima_new[[i]] <- optima_new_df$optima
    int_mat_new[[i]] <- int_mat[new_pool, new_pool]
    
  }
  
  spe_inf <- do.call("rbind", spe_inf)
  
  return(list(spe_inf = spe_inf,
              S_size = S_size,
              optima_new = optima_new,
              int_mat_new = int_mat_new))
  
}
#################################################################################
################################################################################
## equal competition-------------------
################################################################################
mcsim_eq_2 <- function(M = 100, S, dispersal, env_df, optima, int_mat, d, landscape,...){
  # d for looping dispersal gradients
  time_sample <- seq(20, 1200, length.out = 60)
  library(dplyr)
  mcsim <- mcomsimr::simulate_MC(patches = M, species = S[d], dispersal = dispersal[d],
                                 torus = T, env.df =  env_df, 
                                 env_optima = optima[[d]], 
                                 int_mat = int_mat[[d]],
                                 landscape = landscape,...)
  
  res <- mcsim$dynamics.df %>% 
    filter(time %in% time_sample) %>%
    mutate(dispersal = dispersal[d])
  res
}


################################################################################
################################################################################
### 6.1 eq_narrow_stage_2
t1 <- Sys.time()
t1

for (i in 1:30) {
  load(paste("meta_df_eq_narrow_", i, 
             ".Rdata", sep = ""))
  
  change <- MC_change(meta_df_eq_narrow, dispersal_gradient = dispersal_gradient, 
                      optima =  narrow_optima, int_mat = eq_int)
  rm(meta_df_eq_narrow)
  env_df  <-  env_df_stage_2[[i]]
  landscape = landscape_stage_2[[i]]
  
  cl <- makeCluster(35)
  registerDoParallel(cl)
  #(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...)
  {meta_df_eq_narrow_stage2 <-foreach(d=c(1:15),
                                        .combine='rbind') %dopar% {
                                          mcsim_eq_2(M = 100, S = change$S_size, 
                                                       dispersal = dispersal_gradient, 
                                                       env_df = env_df, d,
                                                       landscape = landscape, 
                                                       optima = change$optima_new, 
                                                       int_mat = change$int_mat_new,
                                                       env_niche_breadth = 0.5)
                                        };}
  
  stopCluster(cl)
  meta_df_eq_narrow_stage2 <- meta_df_eq_narrow_stage2 %>% 
    full_join(., change$spe_inf, by = c("species", "dispersal")) %>% mutate(reps = i)
  
  save(meta_df_eq_narrow_stage2, file = paste("stage2/", "meta_df_eq_narrow_stage2_", i, 
                                                ".Rdata", sep = ""))
  print(paste(i,"/30"))
  rm(meta_df_eq_narrow_stage2)
  rm(change)
  #setTxtProgressBar(pb, rep)
}
t2 <- Sys.time()
t2
t2-t1

################################################################################
################################################################################
### 6.2 eq_flat_stage_2
t1 <- Sys.time()
t1

for (i in 1:30) {
  load(paste("meta_df_eq_flat_", i, 
             ".Rdata", sep = ""))
  
  change <- MC_change(meta_df_eq_flat, dispersal_gradient = dispersal_gradient, 
                      optima =  flat_optima, int_mat = eq_int)
  rm(meta_df_eq_flat)
  
  env_df = env_df_stage_2[[i]]
  landscape = landscape_stage_2[[i]]
  
  cl <- makeCluster(35)
  registerDoParallel(cl)
  #(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...)
  {meta_df_eq_flat_stage2 <-foreach(d = 1:15,.combine='rbind') %dopar% {
    mcsim_eq_2(M = 100,S = change$S_size, 
                 dispersal = dispersal_gradient, 
                 env_df = env_df, d,
                 landscape = landscape, 
                 optima = change$optima_new,
                 int_mat = change$int_mat_new,
                 env_niche_breadth = 10)
  };}
  
  stopCluster(cl)
  meta_df_eq_flat_stage2 <- meta_df_eq_flat_stage2 %>% 
    full_join(., change$spe_inf, by = c("species", "dispersal")) %>% 
    mutate(reps = i)
  
  save(meta_df_eq_flat_stage2, file = paste("stage2/", "meta_df_eq_flat_stage2_", i, 
                                              ".Rdata", sep = ""))
  print(paste(i,"/30"))
  rm(meta_df_eq_flat_stage2)
  rm(change)
  #setTxtProgressBar(pb, rep)
}
t2 <- Sys.time()
t2
t2-t1

################################################################################
## stabilizing competition----------------
################################################################################
## function to simulate species and patch change
MC_change <- function(MC_out, dispersal_gradient, optima, int_mat){
  
  MC_t1 <- MC_out %>% 
    filter(time == 1200) %>% # select the last step
    filter(patch %in% 21:100) %>% ## remove the first 20 patches
    filter(N > 0)  # keep species that survive
  #spe_pool <- list()
  S_size <- vector(length = 15)
  spe_inf <- list()
  optima_new <- list()
  int_mat_new <- list()
  
  for (i in 1:15){
    species_left_0 <- unique((MC_t1 %>% 
                                filter(dispersal == dispersal_gradient[i]))$`species`)
    
    species_left <- species_left_0[3:length(species_left_0)] # remove 2 species
    new_pool <- c(species_left, 51:70) # introduce 20 species to meta-community
    S_size[i] <- length(new_pool)
    spe_inf[[i]] <- data.frame(dispersal = dispersal_gradient[i], 
                              species = 1:length(new_pool), species2 = new_pool)
    optima_new_df <- optima %>% 
      filter(species %in% new_pool)
    
    optima_new[[i]] <- optima_new_df$optima
    int_mat_new[[i]] <- int_mat[new_pool, new_pool]
    
  }
  
  spe_inf <- do.call("rbind", spe_inf)
  
  return(list(spe_inf = spe_inf,
              S_size = S_size,
              optima_new = optima_new,
              int_mat_new = int_mat_new))
  
}

################################################################################
mcsim_stab_2 <- function(M = 100, S, dispersal, env_df, optima, int_mat, d, landscape,...){
  # d for loop dispersal gradients, rep for looping  landscape
  time_sample <- seq(20, 1200, length.out = 60)
  library(dplyr)
  mcsim <- mcomsimr::simulate_MC(patches = M, species = S[d], 
                                 dispersal = dispersal[d],
                                 torus = T, env.df =  env_df, 
                                 env_optima = optima[[d]], 
                                 int_mat = int_mat[[d]],
                                 landscape = landscape, ...)
  
  res <- mcsim$dynamics.df %>% 
    filter(time %in% time_sample) %>%
    mutate(dispersal = dispersal[d])
  res
}



################################################################################
################################################################################
### 6.1 stab_narrow_stage_2
t1 <- Sys.time()
t1

for (i in 1:30) {
  
  load(paste("meta_df_stab_narrow_", i, 
             ".Rdata", sep = ""))
  
  change <- MC_change(meta_df_stab_narrow, dispersal_gradient = dispersal_gradient, 
                      optima =  narrow_optima, int_mat = stab_int)
  rm(meta_df_stab_narrow)
  env_df  <-  env_df_stage_2[[i]]
  landscape = landscape_stage_2[[i]]
  
  cl <- makeCluster(35)
  registerDoParallel(cl)
  #(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...)
  {meta_df_stab_narrow_stage2 <-foreach(d=c(1:15),
                                        .combine='rbind') %dopar% {
                                          mcsim_stab_2(M = 100, S = change$S_size, 
                                                       dispersal = dispersal_gradient, 
                                                       env_df = env_df, d,
                                                       landscape = landscape, 
                                                       optima = change$optima_new, 
                                                       int_mat = change$int_mat_new,
                                                       env_niche_breadth = 0.5)
                                        };}
  
  stopCluster(cl)
  meta_df_stab_narrow_stage2 <- meta_df_stab_narrow_stage2 %>% 
    full_join(., change$spe_inf, by = c("species", "dispersal")) %>% mutate(reps = i)
  
  save(meta_df_stab_narrow_stage2, file = paste("stage2/", "meta_df_stab_narrow_stage2_", i, 
                                              ".Rdata", sep = ""))
  print(paste(i,"/30"))
  rm(meta_df_stab_narrow_stage2)
  rm(change)
  #setTxtProgressBar(pb, rep)
}
t2 <- Sys.time()
t2
t2-t1

################################################################################
################################################################################
### 6.2 stab_flat_stage_2
t1 <- Sys.time()
t1

for (i in 1:30) {
  load(paste("meta_df_stab_flat_", i, 
             ".Rdata", sep = ""))
  
  change <- MC_change(meta_df_stab_flat, dispersal_gradient = dispersal_gradient, 
                      optima =  flat_optima, int_mat = stab_int)
  rm(meta_df_stab_flat)
  
  env_df  <-  env_df_stage_2[[i]]
  landscape = landscape_stage_2[[i]]
  
  cl <- makeCluster(35)
  registerDoParallel(cl)
  #(M = 100, S = 50, dispersal, env_df, optima,rep,landscape,...)
  {meta_df_stab_flat_stage2 <-foreach(d=1:15,
                                        .combine='rbind') %dopar% {
                                          mcsim_stab_2(M = 100, S = change$S_size, dispersal = dispersal_gradient, 
                                                       env_df = env_df, d,
                                                       landscape = landscape, 
                                                       optima = change$optima_new, 
                                                       int_mat = change$int_mat_new,
                                                       env_niche_breadth = 10)
                                        };}
  
  stopCluster(cl)
  meta_df_stab_flat_stage2 <- meta_df_stab_flat_stage2 %>% 
    full_join(., change$spe_inf, by = c("species", "dispersal")) %>% mutate(reps = i)
  
  save(meta_df_stab_flat_stage2, file = paste("stage2/", "meta_df_stab_flat_stage2_", i, 
                                                ".Rdata", sep = ""))
  print(paste(i,"/30"))
  rm(meta_df_stab_flat_stage2)
  rm(change)
  #setTxtProgressBar(pb, rep)
}
t2 <- Sys.time()
t2
t2-t1






