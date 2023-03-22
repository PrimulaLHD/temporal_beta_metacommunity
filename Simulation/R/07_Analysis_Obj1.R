##this is the r script to analysis the simulated metacommunity data
## 2022-03-08

################################################################################
## Part 1 functions--------
################################################################################
library(bipartite)
library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(mgcv)
library(patchwork)


### temporal beta diversity of species-patch links
#######new function
beta_temporal_metacomm <- function(m1, m2, method = "Jac") {
  # m1 and m2 are the metacommunity matrix in time 1 and time 2, respectively
  # or species-habitat network in time 1, and in time 2, respectively
  if (all(dim(m1) != dim(m2))) {
    print("Error, m1 and m2 must have equal dimensions!")
  }
  
  if (all(dim(m1) == dim(m2))) {
    # share items
    D <- m1 - m2
    a <- sum(m1[m1 == m2]) #  number of shared items between m1 and m2
    b <- length(D[D == 1]) # number of unique items in m1
    c <- length(D[D == -1]) # number of unique items in m2
    
    
    # part 1, local components
    m1_0 <- bipartite::empty(m1)
    m2_0 <- bipartite::empty(m2)
    
    shared_sp1 <- colnames(m1_0) %in% colnames(m2_0)
    shared_patch1 <- rownames(m1_0) %in% rownames(m2_0)
    
    shared_sp2 <- colnames(m2_0) %in% colnames(m1_0)
    shared_patch2 <- rownames(m2_0) %in% rownames(m1_0)
    
    m1_00 <- (m1_0[shared_patch1, shared_sp1])
    m2_00 <- (m2_0[shared_patch2, shared_sp2])
    
    D0 <- m1_00 - m2_00
    b0 <-
      length(D0[D0 == 1]) # unique numbers of item in m1, local-extinction driven part
    c0 <-
      length(D0[D0 == -1]) # unique numbers of item in m2, local-colonization driven part
    
    # part 2, regional species changes
    
    R_mat_loss <- D[rowSums(m2) > 0, colSums(m2) == 0]
    b_R <-
      length(R_mat_loss[R_mat_loss == 1]) # regional species loss (extinction) driven part
    
    
    R_mat_gain <- D[rowSums(m1) > 0, colSums(m1) == 0]
    c_R <-
      length(R_mat_gain[R_mat_gain == -1]) # regional species gain (colonization) driven part
    
    # part 3, landscape modification
    
    L_mat_loss <- D[rowSums(m2) == 0, colSums(m2) > 0]
    b_L <- length(L_mat_loss[L_mat_loss == 1])
    
    L_mat_gain <- D[rowSums(m1) == 0, colSums(m1) > 0]
    c_L <- length(L_mat_gain[L_mat_gain == -1])
    
    # part 4, both landscape modification and regional species change driven part
    
    RL_mat_loss <- D[rowSums(m2) == 0, colSums(m2) == 0]
    b_RL <- length(RL_mat_loss[RL_mat_loss == 1])
    
    RL_mat_gain <- D[rowSums(m1) == 0, colSums(m1) == 0]
    c_RL <- length(RL_mat_gain[RL_mat_gain == -1])
    
    if (method == "Jac") {
      # Jaccard dissimilarity
      beta_Local <- (b0 + c0) / (a + b + c)
      beta_Regional <- (b_R + c_R) / (a + b + c)
      beta_Landscape <- (b_L + c_L) / (a + b + c)
      beta_RL <- (b_RL + c_RL) / (a + b + c)
      beta_temporal <-
        ((b + c) / (a + b + c)) # total temporal beta diversity of link types
      
      ## further partitioning
      beta_extinction <- (b0 + b_R + b_L + b_RL) / (a + b + c)
      beta_colonization <- (c0 + c_R + c_L + c_RL) / (a + b + c)
      beta_Local_colonization <- c0 / (a + b + c)
      beta_Local_extinction <- b0  / (a + b + c)
      beta_Regional_colonization <- c_R / (a + b + c)
      beta_Regional_extinction <- b_R / (a + b + c)
      beta_Landscape_gain <- c_L / (a + b + c)
      beta_Landscape_loss <- b_L / (a + b + c)
      beta_RL_colonization <- c_RL / (a + b + c)
      beta_RL_extinction <- b_RL / (a + b + c)
      
    }
    
    if (method == "Sor") {
      #Sørensen
      beta_Local <- (b0 + c0) / (2 * a + b + c)
      beta_Regional <- (b_R + c_R) / (2 * a + b + c)
      beta_Landscape <- (b_L + c_L) / (2 * a + b + c)
      beta_RL <- (b_RL + c_RL) / (2 * a + b + c)
      beta_temporal <-
        ((b + c) / (2 * a + b + c)) # total temporal beta diversity of link types
      
      beta_extinction <- (b0 + b_R + b_L + b_RL) / (2 * a + b + c)
      beta_colonization <- (c0 + c_R + c_L + c_RL) / (2 * a + b + c)
      
      ## further partitioning
      beta_extinction <- (b0 + b_R + b_L + b_RL) / (2 * a + b + c)
      beta_colonization <- (c0 + c_R + c_L + c_RL) / (2 * a + b + c)
      beta_Local_colonization <- c0 / (2 * a + b + c)
      beta_Local_extinction <- b0  / (2 * a + b + c)
      beta_Regional_colonization <- c_R / (2 * a + b + c)
      beta_Regional_extinction <- b_R / (2 * a + b + c)
      beta_Landscape_gain <- c_L / (2 * a + b + c)
      beta_Landscape_loss <- b_L / (2 * a + b + c)
      beta_RL_colonization <- c_RL / (2 * a + b + c)
      beta_RL_extinction <- b_RL / (2 * a + b + c)
      
      
    }
    
  }
  
  res <-
    data.frame(
      Component = c(
        "beta_Local",
        "beta_Regional",
        "beta_Landscape",
        "beta_RL",
        "beta_temporal",
        "beta_extinction",
        "beta_colonization",
        "beta_Local_colonization",
        "beta_Local_extinction",
        "beta_Regional_colonization",
        "beta_Regional_extinction",
        "beta_Landscape_gain",
        "beta_Landscape_loss",
        "beta_RL_colonization",
        "beta_RL_extinction"
      ),
      beta = c(
        beta_Local,
        beta_Regional,
        beta_Landscape,
        beta_RL,
        beta_temporal,
        beta_extinction,
        beta_colonization,
        beta_Local_colonization,
        beta_Local_extinction,
        beta_Regional_colonization,
        beta_Regional_extinction,
        beta_Landscape_gain,
        beta_Landscape_loss,
        beta_RL_colonization,
        beta_RL_extinction
      )
    )
  
  return(res)
  
}


### function to make quantitative network to binary network
q2bnet <- function(net){
  net[net > 0] <- 1
  net
}

### function for calculation network metrics for each time step
net_str1 <- function(mcsim_out, dispersal_gradient, niche, competition, i, rep){
  ## i for looping dispersal gradients, and rep for looping landscapes
  library(dplyr)
  webs <- mcsim_out %>% 
    filter(dispersal == dispersal_gradient[i]) %>% 
    bipartite::frame2webs(., varnames = c("patch", "species", "time", "N"))
  
  gamma <- sapply(webs, function(x) length(colnames(x)))
  alpha <- sapply(webs, function(x) mean(rowSums(x>0)))
  links <- sapply(webs, function(x) sum(x>0))
  connectance <- sapply(webs, function(x) bipartite::networklevel(x, index = "connectance"))
  nodf <- sapply(webs, function(x) bipartite::networklevel(x, index = "NODF"))
  #h2 <- sapply(webs, function(x) bipartite::H2fun(x)[1])
  modularity <- sapply(webs, function(x) ifelse(length(colnames(x)) > 1, bipartite::computeModules(x)@likelihood, NA))
  nets_str <- data.frame(links = links, 
                         connectance = connectance,
                         nodf = nodf, 
                         modularity = modularity, 
                         gamma = gamma,
                         alpha = alpha,
                         dispersal = dispersal_gradient[i],
                         reps = rep) 
  nets_str$time <- rownames(nets_str)
  res2 <- reshape2::melt(nets_str, id = c("time", "reps", "dispersal"))
  names(res2)[4:5] <- c("structure", "value")
  res3 <- res2 %>% mutate(niche = niche, competition = competition)
  res3
  
}

### function for calculation temporal beta diversity network metrics for each time series

beta_link1 <- function(mcsim_out, dispersal_gradient, niche, competition, rep){
  library(dplyr)
  beta_temporal <- list()
  beta_ini_t <- list()
  beta_temporal_ini_final <- list()
  div_t <- list()# species diversity in time point t
  for (i in 1:15) {
    webs <- mcsim_out %>% 
      filter(dispersal == dispersal_gradient[i]) %>% 
      bipartite::frame2webs(., varnames = c("patch", "species", "time", "N"), 
                            emptylist = F) 
    betares <- list()
    beta_t <- list() # beta diversity between any time to initial time
    for (j in 1:59) {
      betares[[j]] <- beta_temporal_metacomm(q2bnet(webs[[j]]), q2bnet(webs[[j+1]]))
      beta_t[[j]] <- beta_temporal_metacomm(q2bnet(webs[[1]]), q2bnet(webs[[j+1]]))
    }
    beta_temporal[[i]] <- do.call("rbind", betares) %>% 
      mutate(niche = niche, 
             competition = competition, 
             dispersal = dispersal_gradient[i],
             rep = rep, time = rep(1:59, each = 15))
    
    beta_ini_t[[i]] <- do.call("rbind", beta_t) %>% 
      mutate(niche = niche, 
             competition = competition, 
             dispersal = dispersal_gradient[i],
             rep = rep, time = rep(1:59, each = 15))
    
    beta_temporal_ini_final[[i]] <- beta_temporal_metacomm(q2bnet(webs[[1]]), 
                                                           q2bnet(webs[[60]])) %>% 
      mutate(niche = niche, 
             competition = competition, 
             dispersal = dispersal_gradient[i],
             rep = rep)
    div <- list()
    for (m in 1:60) {
      div[[m]] <- data.frame(gamma_richness = length(colnames(empty(webs[[m]]))),
                             alpha_richness = mean(rowSums((empty(webs[[m]])) > 0)),
                             time = m)
    }
    
    div_t[[i]] <- do.call("rbind", div) %>% 
      mutate(niche = niche, 
             competition = competition, 
             dispersal = dispersal_gradient[i],
             rep = rep)
    
    
  }
  
  beta_temporal <- do.call("rbind", beta_temporal)
  beta_temporal_ini_t <- do.call("rbind", beta_ini_t)
  beta_temporal_ini_final <- do.call("rbind", beta_temporal_ini_final)
  div_t <- do.call("rbind", div_t)
  
  return(list(beta_temporal,
              beta_temporal_ini_t,
              beta_temporal_ini_final,
              div_t))
  
}

## colors
mycol <- RColorBrewer::brewer.pal(7,"Dark2")


################################################################################
### part 2 data analysis for object1---------
################################################################################
library(bipartite)
library(mcomsimr)
library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(mgcv)
library(patchwork)
library(cowplot)

setwd("/Users/haidong/dropbox/GitHub/temporal_beta_metacommunity")

##########################################################################
# Object 1 the relationship between processes and emerged network----
##########################################################################
## 1. calculations
### 1.1 equal competition scenarios-------------

#### 1.1.1 equal competition with narrow niche-------

## network structure-----
t1 <- Sys.time()
netstr_eq_narrow <- list()
for (rep in 1:30) {
  load(paste("Data/data2/meta_df_eq_narrow_", rep,".Rdata",sep = ""))
  cl <- makeCluster(35)
  registerDoParallel(cl)
  {netstr_eq_narrow[[rep]] <-foreach(i = 1:15, .combine='rbind') %dopar% {
    net_str1(meta_df_eq_narrow, dispersal_gradient = dispersal_gradient, 
             niche = "narrow", competition = "equal_competition", i, rep = rep)
  };}
  
  stopCluster(cl)
  rm(meta_df_eq_narrow)
  print(paste("finished ", rep,"/30 of network structure calculation runs!", sep = ""))
}
netstr_eq_narrow2 <- do.call("rbind", netstr_eq_narrow)

t2 <- Sys.time()
t2-t1

save(netstr_eq_narrow2, file = "Res/R1/netstr_eq_narrow2.Rata")
load("Res/R1/netstr_eq_narrow2.Rata")

################################################################################
## beta diversity--------
beta_links_eq_narrow <- list()
beta_links_eq_narrow_ini_t <- list()
beta_links_eq_narrow_ini_finnal <- list()
div_eq_narrow <- list()
t1 <- Sys.time()
for (rep in 1:30) {
  load(paste("Data/data2/meta_df_eq_narrow_", rep,".Rdata",sep = ""))
  tep <- beta_link1(meta_df_eq_narrow, 
                    dispersal_gradient = dispersal_gradient, 
                    niche = "narrow", competition="equal_competition", rep = rep)
  beta_links_eq_narrow[[rep]] <- tep[[1]]
  beta_links_eq_narrow_ini_t[[rep]] <- tep[[2]]
  beta_links_eq_narrow_ini_finnal[[rep]] <- tep[[3]]
  div_eq_narrow[[rep]] <- tep[[4]]
  
  rm(meta_df_eq_narrow)
  print(paste("finished ", rep,"/30 of network beta links calculation runs!", sep = ""))
}

t2 <- Sys.time()
t2-t1

beta_links_eq_narrow  <- do.call("rbind", beta_links_eq_narrow)


#### 1.1.2 equal competition with flat niche -----------------------------------

## network structure----
t1 <- Sys.time()
netstr_eq_flat <- list()
for (rep in 1:30) {
  load(paste("Data/data2/meta_df_eq_flat_", rep,".Rdata",sep = ""))
  cl <- makeCluster(35)
  registerDoParallel(cl)
  {netstr_eq_flat[[rep]] <- foreach(i = 1:15, .combine='rbind') %dopar% {
    net_str1(meta_df_eq_flat, dispersal_gradient = dispersal_gradient, 
             niche = "flat", competition = "equal_competition", i, rep = rep)
  };}
  
  stopCluster(cl)
  rm(meta_df_eq_flat)
  print(paste("finished ", rep,"/30 of network structure calculation runs!", sep = ""))
}
netstr_eq_flat2_2 <- do.call("rbind", netstr_eq_flat)

t2 <- Sys.time()
t2-t1

save(netstr_eq_flat2_2, file = "Res/R1/netstr_eq_flat2.Rata")

load("Res/R1/netstr_eq_flat2.Rata")


################################################################################
## beta diversity--------
beta_links_eq_flat <- list()
beta_links_eq_flat_ini_t <- list()
beta_links_eq_flat_ini_finnal <- list()
div_eq_flat <- list()
t1 <- Sys.time()
t1
for (rep in 1:30) {
  load(paste("Data/data2/meta_df_eq_flat_", rep,".Rdata",sep = ""))
  tep <- beta_link1(meta_df_eq_flat, 
                    dispersal_gradient = dispersal_gradient, 
                    niche = "flat", competition="equal_competition", rep = rep)
  beta_links_eq_flat[[rep]] <- tep[[1]]
  beta_links_eq_flat_ini_t[[rep]] <- tep[[2]]
  beta_links_eq_flat_ini_finnal[[rep]] <- tep[[3]]
  div_eq_flat[[rep]] <- tep[[4]]
  
  rm(meta_df_eq_flat)
  print(paste("finished ", rep,"/30 of network beta links calculation runs!", sep = ""))
}

t2 <- Sys.time()
t2-t1

beta_links_eq_flat  <- do.call("rbind", beta_links_eq_flat)





## 1.2 stabilizing competition scenarios----------------------------------------
### 1.2.1 stabilizing competition with narrow niche---------------------------

### network structure----
t1 <- Sys.time()
netstr_stab_narrow <- list()
for (rep in 1:30) {
  load(paste("Data/data2/meta_df_stab_narrow_", rep,".Rdata",sep = ""))
  
  cl <- makeCluster(35)
  registerDoParallel(cl)
  {netstr_stab_narrow[[rep]] <-foreach(i = 1:15, .combine='rbind') %dopar% {
    net_str1(meta_df_stab_narrow, dispersal_gradient = dispersal_gradient, 
             niche = "narrow", competition = "stable_competition", i, rep = rep)
  };}
  
  stopCluster(cl)
  rm(meta_df_stab_narrow)
  print(paste("finished ", rep,"/30 of network structure calculation runs!", sep = ""))
}
netstr_stab_narrow2 <- do.call("rbind", netstr_stab_narrow)

t2 <- Sys.time()
t2-t1

save(netstr_stab_narrow2, file = "Res/R1/netstr_stab_narrow2.Rata")
load("Res/R1/netstr_stab_narrow2.Rata")


### beta diversity--------
beta_links_stab_narrow <- list()
beta_links_stab_narrow_ini_t <- list()
beta_links_stab_narrow_ini_finnal <- list()
div_stab_narrow <- list()
t1 <- Sys.time()
for (rep in 1:30) {
  load(paste("Data/data2/meta_df_stab_narrow_", rep,".Rdata",sep = ""))
  tep <- beta_link1(meta_df_stab_narrow, 
                    dispersal_gradient = dispersal_gradient, 
                    niche = "narrow", competition = "stabilizing_competition", rep = rep)
  beta_links_stab_narrow[[rep]] <- tep[[1]]
  beta_links_stab_narrow_ini_t[[rep]] <- tep[[2]]
  beta_links_stab_narrow_ini_finnal[[rep]] <- tep[[3]]
  div_stab_narrow[[rep]] <- tep[[4]]
  
  rm(meta_df_stab_narrow)
  print(paste("finished ", rep,"/30 of network beta links calculation runs!", sep = ""))
}

t2 <- Sys.time()
t2-t1

beta_links_stab_narrow  <- do.call("rbind", beta_links_stab_narrow)

#################
### 1.2.2 stabilizing competition with flat niche-------

### network structure------
t1 <- Sys.time()
netstr_stab_flat <- list()
for (rep in 1:30) {
  load(paste("Data/data2/meta_df_stab_flat_", rep,".Rdata",sep = ""))
  
  cl <- makeCluster(35)
  registerDoParallel(cl)
  {netstr_stab_flat[[rep]] <-foreach(i = 1:15, .combine='rbind') %dopar% {
    net_str1(meta_df_stab_flat, dispersal_gradient = dispersal_gradient, 
             niche = "flat", competition = "stable_competition", i, rep = rep)
  };}
  
  stopCluster(cl)
  rm(meta_df_stab_flat)
  print(paste("finished ", rep,"/30 of network structure calculation runs!", sep = ""))
}
netstr_stab_flat2 <- do.call("rbind", netstr_stab_flat)

t2 <- Sys.time()
t2-t1

save(netstr_stab_flat2, file = "Res/R1/netstr_stab_flat2.Rata")
load("Res/R1/netstr_stab_flat2.Rata")

# load("Res/R1/netstr_stab_flat2_2.Rata")
# 
# netstr_stab_flat2_2 <- netstr_stab_flat2
# rm(netstr_stab_flat2)
# 
# netstr_stab_flat2 <- rbind(netstr_stab_flat2_2, netstr_stab_flat2)
# rm(netstr_stab_flat2_2)


### beta diversity

### beta diversity--------
beta_links_stab_flat <- list()
beta_links_stab_flat_ini_t <- list()
beta_links_stab_flat_ini_finnal <- list()

t1 <- Sys.time()
for (rep in 1:30) {
  load(paste("Data/data2/meta_df_stab_flat_", rep,".Rdata",sep = ""))
  tep <- beta_link1(meta_df_stab_flat, 
                    dispersal_gradient = dispersal_gradient, 
                    niche = "flat", competition = "stabilizing_competition", rep = rep)
  beta_links_stab_flat[[rep]] <- tep[[1]]
  beta_links_stab_flat_ini_t[[rep]] <- tep[[2]]
  beta_links_stab_flat_ini_finnal[[rep]] <- tep[[3]]
  
  
  rm(meta_df_stab_flat)
  print(paste("finished ", rep,"/30 of network beta links calculation runs!", sep = ""))
}

t2 <- Sys.time()
t2-t1

beta_links_stab_flat  <- do.call("rbind", beta_links_stab_flat)



################################################################################
## 2.make figures---------------------------------------------------------------
## network structures
net_Res <- do.call("rbind", list(
  na.omit(netstr_eq_narrow2) %>% 
    group_by(reps, dispersal, structure, niche, competition) %>%
    summarise(value = mean(value)) %>%
    group_by(dispersal, structure, niche, competition) %>% 
    summarise(value_m = median(value),
              value_h = quantile(value, 0.75),
              value_l = quantile(value, 0.25)) %>% 
    mutate(competition = "equal competition",
           niche = "narrow"),
  na.omit(netstr_eq_flat2) %>% 
    group_by(reps, dispersal, structure, niche, competition) %>%
    summarise(value = mean(value)) %>%
    group_by(dispersal, structure, niche, competition) %>% 
    summarise(value_m = median(value),
              value_h = quantile(value, 0.75),
              value_l = quantile(value, 0.25)) %>% 
    mutate(competition = "equal competition",
           niche = "flat"),
  na.omit(netstr_stab_narrow2) %>% 
    group_by(reps, dispersal, structure, niche, competition) %>%
    summarise(value = mean(value)) %>%
    group_by(dispersal, structure, niche, competition) %>% 
    summarise(value_m = median(value),
              value_h = quantile(value, 0.75),
              value_l = quantile(value, 0.25)) %>% 
    mutate(competition = "stabilizing competition",
           niche = "narrow"),
  na.omit(netstr_stab_flat2) %>% 
    group_by(reps, dispersal, structure, niche, competition) %>%
    summarise(value = mean(value)) %>%
    group_by(dispersal, structure, niche, competition) %>% 
    summarise(value_m = median(value),
              value_h = quantile(value, 0.75),
              value_l = quantile(value, 0.25)) %>% 
    mutate(competition = "stabilizing competition",
           niche = "flat")
))

net_Res$competition <- factor(net_Res$competition, levels = c("equal competition", 
                                                              "stabilizing competition"),
                              labels = c("equal competition", 
                                         "stabilizing competition"))
fig3_a <- ggplot(net_Res,
                 aes(x = dispersal, y = value_m, fill = factor(structure,
                                                               levels = c("links", "connectance",
                                                                          "nodf", "modularity",
                                                                          "gamma", "alpha"),
                                                               labels = c("Number of links", 
                                                                          "Connectance",
                                                                          "Nestedness, NODF", 
                                                                          "Modularity, M",
                                                                          "gamma-richness", 
                                                                          "alpha-richness")),
                     color = factor(structure,
                                    levels = c("links", "connectance",
                                               "nodf", "modularity",
                                               "gamma", "alpha"),
                                    labels = c("Number of links", 
                                               "Connectance",
                                               "Nestedness, NODF", "Modularity, M",
                                               "gamma-richness", 
                                               "alpha-richness"))))+
  geom_ribbon(aes(ymin = value_l, ymax = value_h), color = NA, alpha = 0.25) +
  geom_line(size = 1) + 
  facet_grid(competition ~  factor(niche, levels = c("narrow", "flat"), 
                                   labels = c("Narrow abiotic niche (\u03C3 = 0.5)",
                                              "Wide abiotic niche (\u03C3 = 10)")))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_color_manual(values = mycol[1:6]) + 
  scale_fill_manual(values = mycol[1:6]) + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Species richness and network structure") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.position = "bottom",
        axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 12,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())


ggsave("./Figs/R1/Figure_network_str_facet_3.png", width = 9*0.8, height = 12*0.60)
ggsave("./Figs/R1/Figure_network_str_facet_3.pdf", width = 9*0.8, height = 12*0.60)


# ggplot(net_Res %>% filter(structure == 'gamma'| structure=='alpha'),
#        aes(x = dispersal, y = value_m, fill = factor(structure),
#            color = factor(structure)))+
#   geom_ribbon(aes(ymin = value_l, ymax = value_h), color = NA, alpha = 0.25) +
#   geom_line(size = 1) + 
#   facet_grid(competition ~  factor(niche, levels = c("narrow", "flat"), 
#                                    labels = c("narrow abiotic niche (\u03C3 = 0.5)",
#                                               "wide abiotic niche (\u03C3 = 10)")))+
#   scale_y_log10()+
#   scale_color_manual(values = mycol[1:6]) +
#   scale_fill_manual(values = mycol[1:6]) + 
#   scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
#   xlab(expression(paste("dispersal (",italic(a[i]),")",sep = ""))) + 
#   ylab("species richness") +
#   theme_bw()+
#   theme(panel.grid = element_blank())+
#   theme(legend.position = "bottom",
#         axis.title.y = element_text(size = 12,family = "serif"),
#         axis.text.y = element_text(size = 12, family = "serif"),
#         axis.text.x = element_text(size = 12, family = "serif"),
#         axis.title.x = element_text(size = 12,family = "serif"),
#         strip.text = element_text(size = 12,family = "serif"),
#         legend.text = element_text(size = 12,family = "serif"),
#         legend.title = element_blank())
# 


## beta diversity plot-
beta_eq_and_stab <- do.call("rbind", list(beta_links_eq_narrow,
                                          beta_links_eq_flat,
                                          beta_links_stab_narrow,
                                          beta_links_stab_flat))

beta_eq_and_stab$competition <- factor(beta_eq_and_stab$competition, 
                                       levels = c("equal_competition", 
                                                  "stabilizing_competition"),
                                       labels = c("equal competition", 
                                                  "stabilizing competition"))


fig3_b <- ggplot(beta_eq_and_stab %>% 
                   filter(Component == "beta_temporal"|
                            Component == "beta_colonization"| 
                            Component == "beta_extinction") %>% 
                   group_by(rep, dispersal, Component, niche, competition) %>% 
                   summarise(value = mean(beta)) %>% 
                   group_by(dispersal, Component, niche, competition) %>% 
                   summarise(value_m = median(value),
                             value_h = quantile(value, 0.75),
                             value_l = quantile(value, 0.25)),
                 aes(x = dispersal, y = value_m, 
                     fill = factor(Component, levels = c("beta_temporal", "beta_colonization", "beta_extinction"),
                                   labels = c("βtemporal", "βcolonization","βextinction")), 
                     color = factor(Component, levels = c("beta_temporal", "beta_colonization", "beta_extinction"),
                                    labels = c("βtemporal", "βcolonization","βextinction")))) +
  geom_ribbon(aes(ymin = value_l, ymax = value_h), color = NA, alpha = 0.25) +
  geom_line(size = 1) + 
  scale_color_manual(values = mycol[1:3]) + 
  scale_fill_manual(values = mycol[1:3]) + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Temporal beta diversity of links") + 
  facet_grid(competition ~  factor(niche, levels = c("narrow", "flat"), 
                                   labels = c("Narrow abiotic niche (\u03C3 = 0.5)",
                                              "Wide abiotic niche (\u03C3 = 10)")))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.position = "bottom",
    axis.title.y=element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 12,family = "serif"),
        legend.text = element_text(size = 12, family = "serif", face = "italic"),
        legend.title = element_blank())

# ggsave("./Figs/Figure_beta_dispersal.png", width = 9*0.8, height = 12*0.60)
# ggsave("./Figs/Figure_beta_dispersal.pdf", width = 9*0.8, height = 12*0.60)

library(patchwork)
fig3_a + fig3_b + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(', tag_sep = '', 
                  tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 12,family = "serif", hjust = 0, vjust = 0))

ggsave("./Figs/R1/Figure_3.png", width = 9*1.2, height = 12*0.55)
ggsave("./Figs/R1/Figure_3.pdf", width = 9*1.2, height = 12*0.55)

################################################################################
## 3. statistics analysis----------
## to test how local and regional diversity effects on emerged network structure
library(data.table)
library(mgcv)
library(visreg)
## equal competition------------
netstr_eq_narrow_lm <- data.table::dcast(netstr_eq_narrow2 %>% 
                                           group_by(reps, dispersal, structure, niche, competition) %>% 
                                           summarise(value = mean(value,na.rm = T)), 
                                         reps + dispersal +niche + competition~structure) 

netstr_eq_narrow_lm$reps <- as.factor(netstr_eq_narrow_lm$reps)
gam_mod_connectance_eq_narrow <- gam(connectance ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                     data = netstr_eq_narrow_lm)
summary(gam_mod_connectance_eq_narrow)

op <- par(mfrow = c(1,2))
visreg(gam_mod_connectance_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_connectance_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)



gam_mod_nodf_eq_narrow <- gam(nodf ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                              data = netstr_eq_narrow_lm)
summary(gam_mod_nodf_eq_narrow)

op <- par(mfrow = c(1,2))
visreg(gam_mod_nodf_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_nodf_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)


gam_mod_modularity_eq_narrow <- gam(modularity ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                    data = netstr_eq_narrow_lm)
summary(gam_mod_modularity_eq_narrow)

op <- par(mfrow = c(1,2))
visreg(gam_mod_modularity_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_modularity_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)


netstr_eq_flat_lm <- data.table::dcast(netstr_eq_flat2 %>% 
                                         group_by(reps, dispersal, structure, niche, competition) %>% 
                                         summarise(value = mean(value,na.rm = T)), 
                                       reps + dispersal +niche + competition~structure) 

netstr_eq_flat_lm$reps <- as.factor(netstr_eq_flat_lm$reps)
gam_mod_connectance_eq_flat <- gam(connectance ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                   data = netstr_eq_flat_lm)
summary(gam_mod_connectance_eq_flat)

op <- par(mfrow = c(1,2))
visreg(gam_mod_connectance_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_connectance_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)



gam_mod_nodf_eq_flat <- gam(nodf ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                            data = netstr_eq_flat_lm)
summary(gam_mod_nodf_eq_flat)

op <- par(mfrow = c(1,2))
visreg(gam_mod_nodf_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_nodf_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)


gam_mod_modularity_eq_flat <- gam(modularity ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                  data = netstr_eq_flat_lm)
summary(gam_mod_modularity_eq_flat)

op <- par(mfrow = c(1,2))
visreg(gam_mod_modularity_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_modularity_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)


## equal competition, net~richness-----
op <- par(mfrow = c(3,4))
visreg(gam_mod_connectance_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_connectance_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

visreg(gam_mod_connectance_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_connectance_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

visreg(gam_mod_nodf_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_nodf_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

visreg(gam_mod_nodf_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_nodf_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
visreg(gam_mod_modularity_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_modularity_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

visreg(gam_mod_modularity_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_modularity_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

par(op)


## stable competition------------

netstr_stab_narrow_lm <- data.table::dcast(netstr_stab_narrow2 %>% 
                                             group_by(reps, dispersal, structure, niche, competition) %>% 
                                             summarise(value = mean(value,na.rm = T)), 
                                           reps + dispersal +niche + competition~structure) 

netstr_stab_narrow_lm$reps <- as.factor(netstr_stab_narrow_lm$reps)
gam_mod_connectance_stab_narrow <- gam(connectance ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                       data = netstr_stab_narrow_lm)
summary(gam_mod_connectance_stab_narrow)

op <- par(mfrow = c(1,2))
visreg(gam_mod_connectance_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_connectance_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)



gam_mod_nodf_stab_narrow <- gam(nodf ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                data = netstr_stab_narrow_lm)
summary(gam_mod_nodf_stab_narrow)

op <- par(mfrow = c(1,2))
visreg(gam_mod_nodf_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_nodf_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)


gam_mod_modularity_stab_narrow <- gam(modularity ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                      data = netstr_stab_narrow_lm)
summary(gam_mod_modularity_stab_narrow)

op <- par(mfrow = c(1,2))
visreg(gam_mod_modularity_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_modularity_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)


netstr_stab_flat_lm <- data.table::dcast(netstr_stab_flat2 %>% 
                                           group_by(reps, dispersal, structure, niche, competition) %>% 
                                           summarise(value = mean(value,na.rm = T)), 
                                         reps + dispersal +niche + competition~structure) 

netstr_stab_flat_lm$reps <- as.factor(netstr_stab_flat_lm$reps)
gam_mod_connectance_stab_flat <- gam(connectance ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                     data = netstr_stab_flat_lm)
summary(gam_mod_connectance_stab_flat)

op <- par(mfrow = c(1,2))
visreg(gam_mod_connectance_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_connectance_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)



gam_mod_nodf_stab_flat <- gam(nodf ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                              data = netstr_stab_flat_lm)
summary(gam_mod_nodf_stab_flat)

op <- par(mfrow = c(1,2))
visreg(gam_mod_nodf_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_nodf_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)


gam_mod_modularity_stab_flat <- gam(modularity ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
                                    data = netstr_stab_flat_lm)
summary(gam_mod_modularity_stab_flat)

op <- par(mfrow = c(1,2))
visreg(gam_mod_modularity_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_modularity_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
par(op)


## stabilizing competition, net~richness----
op <- par(mfrow = c(3,4))
visreg(gam_mod_connectance_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_connectance_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

visreg(gam_mod_connectance_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_connectance_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

visreg(gam_mod_nodf_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_nodf_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

visreg(gam_mod_nodf_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_nodf_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")
visreg(gam_mod_modularity_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_modularity_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

visreg(gam_mod_modularity_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness")

visreg(gam_mod_modularity_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness")

par(op)


#####################################################################################
### beta-diversity----------

## beta_eq_narrow----
beta_eq_narrow_lm <- data.table::dcast(netstr_eq_narrow2 %>% 
                                         group_by(reps, dispersal, structure, niche, competition) %>% 
                                         summarise(value = mean(value,na.rm = T)), 
                                       reps + dispersal +niche + competition~structure) %>%
  full_join(.,beta_links_eq_narrow %>% 
              filter(Component == "beta_temporal"|
                       Component == "beta_colonization"| 
                       Component == "beta_extinction") %>% 
              group_by(rep, dispersal, Component, niche, competition) %>% 
              summarise(beta = mean(beta)) %>% 
              rename(reps = rep), 
            by = c("reps", "dispersal", "niche", "competition")) 

beta_eq_narrow_lm$reps <- as.factor(beta_eq_narrow_lm$reps)

gam_mod_temporal_eq_narrow <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_eq_narrow_lm,
                    Component == "beta_temporal"))

summary(gam_mod_temporal_eq_narrow)

vis.gam(gam_mod_temporal_eq_narrow)
visreg(gam_mod_temporal_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_temporal"))


gam_mod_colonization_eq_narrow <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_eq_narrow_lm,
                    Component == "beta_colonization"))

summary(gam_mod_colonization_eq_narrow)

vis.gam(gam_mod_colonization_eq_narrow)
visreg(gam_mod_colonization_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "colonization component",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points = list(cex = 1), xlab = "gamma-richness",
       ylab = "colonization component",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_colonization"))

gam_mod_extinction_eq_narrow <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_eq_narrow_lm,
                    Component == "beta_extinction"))

summary(gam_mod_extinction_eq_narrow)

vis.gam(gam_mod_extinction_eq_narrow)
visreg(gam_mod_extinction_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "extinction component",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_extinction"))

visreg(gam_mod_extinction_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "extinction component",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_extinction"))

## beta_eq_flat----
beta_eq_flat_lm <- data.table::dcast(netstr_eq_flat2 %>% 
                                       group_by(reps, dispersal, structure, niche, competition) %>% 
                                       summarise(value = mean(value,na.rm = T)), 
                                     reps + dispersal +niche + competition~structure) %>%
  full_join(beta_links_eq_flat %>% 
              filter(Component == "beta_temporal"|
                       Component == "beta_colonization"| 
                       Component == "beta_extinction") %>% 
              group_by(rep, dispersal, Component, niche, competition) %>% 
              summarise(beta = mean(beta)) %>% 
              rename(reps = rep), 
            by = c("reps", "dispersal", "niche", "competition")) 
beta_eq_flat_lm$reps <- as.factor(beta_eq_flat_lm$reps)

gam_mod_temporal_eq_flat <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_eq_flat_lm,
                    Component == "beta_temporal"))

summary(gam_mod_temporal_eq_flat)

vis.gam(gam_mod_temporal_eq_flat)
visreg(gam_mod_temporal_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_temporal"))


gam_mod_colonization_eq_flat <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_eq_flat_lm,
                    Component == "beta_colonization"))

summary(gam_mod_colonization_eq_flat)

vis.gam(gam_mod_colonization_eq_flat)

visreg(gam_mod_colonization_eq_flat,"alpha", scale = "response", 
       partial = TRUE,
       points = list(cex=1), xlab = "alpha-richness",
       ylab = "colonization component",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_eq_flat,"gamma", scale = "response", 
       partial = TRUE,
       points = list(cex=1), xlab = "gamma-richness",
       ylab = "colonization component",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_colonization"))

gam_mod_extinction_eq_flat <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_eq_flat_lm,
                    Component == "beta_extinction"))

summary(gam_mod_extinction_eq_flat)

vis.gam(gam_mod_extinction_eq_flat)
visreg(gam_mod_extinction_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "extinction component",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_extinction"))

visreg(gam_mod_extinction_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "extinction component",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_extinction"))


## beta_stab_narrow----
beta_stab_narrow_lm <- data.table::dcast(netstr_stab_narrow2 %>% 
                                           group_by(reps, dispersal, structure, niche, competition) %>% 
                                           summarise(value = mean(value,na.rm = T)), 
                                         reps + dispersal +niche + competition~structure) %>%
  full_join(.,beta_links_stab_narrow %>% 
              filter(Component == "beta_temporal"|
                       Component == "beta_colonization"| 
                       Component == "beta_extinction") %>% 
              group_by(rep, dispersal, Component, niche, competition) %>% 
              summarise(beta = mean(beta)) %>% 
              mutate(competition = "stable_competition") %>% 
              rename(reps = rep), 
            by = c("reps", "dispersal", "niche", "competition")) 

beta_stab_narrow_lm$reps <- as.factor(beta_stab_narrow_lm$reps)

gam_mod_temporal_stab_narrow <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_stab_narrow_lm,
                    Component == "beta_temporal"))

summary(gam_mod_temporal_stab_narrow)

vis.gam(gam_mod_temporal_stab_narrow)
visreg(gam_mod_temporal_stab_narrow,"alpha", scale = "response", 
       partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_temporal"))


gam_mod_colonization_stab_narrow <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_stab_narrow_lm,
                    Component == "beta_colonization"))

summary(gam_mod_colonization_stab_narrow)

vis.gam(gam_mod_colonization_stab_narrow)
visreg(gam_mod_colonization_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "colonization component",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "colonization component",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_colonization"))

gam_mod_extinction_stab_narrow <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_stab_narrow_lm,
                    Component == "beta_extinction"))

summary(gam_mod_extinction_stab_narrow)

vis.gam(gam_mod_extinction_stab_narrow)
visreg(gam_mod_extinction_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "extinction component",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_extinction"))

visreg(gam_mod_extinction_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "extinction component",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_extinction"))

## beta_stab_flat----
beta_stab_flat_lm <- data.table::dcast(netstr_stab_flat2 %>% 
                                         group_by(reps, dispersal, structure, niche, competition) %>% 
                                         summarise(value = mean(value,na.rm = T)), 
                                       reps + dispersal +niche + competition~structure) %>%
  full_join(beta_links_stab_flat %>% 
              filter(Component == "beta_temporal"|
                       Component == "beta_colonization"| 
                       Component == "beta_extinction") %>% 
              group_by(rep, dispersal, Component, niche, competition) %>% 
              summarise(beta = mean(beta)) %>% 
              mutate(competition = "stable_competition") %>% 
              rename(reps = rep), 
            by = c("reps", "dispersal", "niche", "competition")) 
beta_stab_flat_lm$reps <- as.factor(beta_stab_flat_lm$reps)

gam_mod_temporal_stab_flat <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_stab_flat_lm,
                    Component == "beta_temporal"))

summary(gam_mod_temporal_stab_flat)

vis.gam(gam_mod_temporal_stab_flat)
visreg(gam_mod_temporal_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_temporal"))


gam_mod_colonization_stab_flat <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_stab_flat_lm,
                    Component == "beta_colonization"))

summary(gam_mod_colonization_stab_flat)

vis.gam(gam_mod_colonization_stab_flat)
visreg(gam_mod_colonization_stab_flat,"alpha", scale = "response", partial = TRUE,
       points = list(cex=1), xlab = "alpha-richness",
       ylab = "colonization component",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_stab_flat,"gamma", scale = "response", partial = TRUE,
       points = list(cex=1), xlab = "gamma-richness",
       ylab = "colonization component",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_colonization"))

gam_mod_extinction_stab_flat <- 
  gam(beta ~ s(alpha) + s(gamma) + s(reps, bs = "re"), 
      data = subset(beta_stab_flat_lm,
                    Component == "beta_extinction"))

summary(gam_mod_extinction_stab_flat)

vis.gam(gam_mod_extinction_stab_flat)
visreg(gam_mod_extinction_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "extinction component",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_extinction"))

visreg(gam_mod_extinction_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "extinction component",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_extinction"))


#### figure_eq_beta_vs_richness-------
op <- par(mfrow = c(3, 4))
visreg(gam_mod_temporal_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_eq_flat,"alpha", scale = "response", partial = TRUE,
       points = list(cex = 1), xlab = "alpha-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_eq_flat,"gamma", scale = "response", partial = TRUE,
       points = list(cex=1), xlab = "gamma-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_colonization_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "colonization component",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_eq_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "colonization component",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_eq_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "colonization component",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "colonization component",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_extinction_eq_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "extinction component",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_extinction"))

visreg(gam_mod_extinction_eq_narrow,"gamma", scale = "response", 
       partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "extinction component",
       data = subset(beta_eq_narrow_lm,
                     Component == "beta_extinction"))

visreg(gam_mod_extinction_eq_flat,"alpha", scale = "response", 
       partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "extinction component",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_extinction"))



visreg(gam_mod_extinction_eq_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "extinction component",
       data = subset(beta_eq_flat_lm,
                     Component == "beta_extinction"))

par(op)


#### figure_stab_beta_vs_richness-------
op <- par(mfrow = c(3, 4))

visreg(gam_mod_temporal_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points = list(cex=1), xlab = "alpha-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_temporal_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "temporal beta diversity of links",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_temporal"))

visreg(gam_mod_colonization_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "colonization component",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "colonization component",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_colonization"))
visreg(gam_mod_colonization_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "colonization component",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_colonization"))

visreg(gam_mod_colonization_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "colonization component",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_colonization"))


visreg(gam_mod_extinction_stab_narrow,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "extinction component",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_extinction"))

visreg(gam_mod_extinction_stab_narrow,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "extinction component",
       data = subset(beta_stab_narrow_lm,
                     Component == "beta_extinction"))


visreg(gam_mod_extinction_stab_flat,"alpha", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "alpha-richness",
       ylab = "extinction component",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_extinction"))

visreg(gam_mod_extinction_stab_flat,"gamma", scale = "response", partial = TRUE,
       points=list(cex=1), xlab = "gamma-richness",
       ylab = "extinction component",
       data = subset(beta_stab_flat_lm,
                     Component == "beta_extinction"))
par(op)

