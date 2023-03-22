################################################################################
# part 3 object 2 temporal changes  between stag1 and 2 -----
################################################################################

## functions to calculate the network structure and temporal beta diversity
#### function eq_narrow_stage
eq_narrow_stage <- function(dir1, dir2, dispersal_gradient) {
  beta_res <- list()
  net_metrics_res <- list()
  
  for (i in 1:30) {
    load(paste(dir1, i, ".Rdata", sep = ""))
    stage1 <- meta_df_eq_narrow %>% 
      filter(time == 1200 & N > 0) %>% 
      select(N, patch, species, time, dispersal) %>% 
      mutate(stage = "t1")
    rm(meta_df_eq_narrow)
    
    load(paste(dir2, i, ".Rdata", sep = ""))
    stage2 <- meta_df_eq_narrow_stage2 %>% 
      filter(time == 1200 & N > 0) %>% 
      select(N, patch, species2, time, dispersal) %>% 
      mutate(stage = "t2", patch = patch + 20) %>% 
      rename(species = species2)
    rm(meta_df_eq_narrow_stage2)
    
    df <- rbind(stage1, stage2)
    beta <- list()
    net_metrics <- list()
    for (m in 1:length(dispersal_gradient)) {
      df_tep <- df %>% filter(dispersal == dispersal_gradient[m])
      webs <- bipartite::frame2webs(df_tep, varnames = c("patch", "species", "stage", "N"),
                                    emptylist = FALSE)
      beta[[m]] <- beta_temporal_metacomm(q2bnet(webs[[1]]), q2bnet(webs[[2]])) %>% 
        mutate(competition = "equal", 
               niche = "narrow", 
               dispersal = dispersal_gradient[m])
      
      gamma <- sapply(webs, function(x) length(colnames(empty(x))))
      alpha <- sapply(webs, function(x) mean(rowSums(empty(x) > 0)))
      links <- sapply(webs, function(x) sum(empty(x) > 0))
      connectance <- sapply(webs, function(x) bipartite::networklevel(empty(x), index = "connectance"))
      nodf <- sapply(webs, function(x) bipartite::networklevel(empty(x), index = "NODF"))
      modularity <- sapply(webs, function(x) ifelse(length(colnames(empty(x))) > 1, 
                                                    bipartite::computeModules(empty(x))@likelihood, NA))
      
      nets_str <- data.frame(links = links, 
                             connectance = connectance,
                             nodf = nodf, 
                             modularity = modularity, 
                             gamma = gamma,
                             alpha = alpha) 
      #print(nets_str)
      nets_str$stage <- rownames(nets_str)
      res2 <- reshape2::melt(nets_str, id = "stage")
      names(res2)[2:3] <- c("structure", "value")
      net_metrics[[m]] <- res2 %>% 
        mutate(competition = "equal", 
               niche = "narrow", 
               dispersal = dispersal_gradient[m])
      
    }
    
    beta_res[[i]] <- do.call("rbind", beta) %>% mutate(rep = i)
    net_metrics_res[[i]] <- do.call("rbind", net_metrics) %>% mutate(rep = i)
  }
  
  return(list(beta = do.call("rbind", beta_res), 
              net_metrics = do.call("rbind", net_metrics_res)))
}


#### function eq_flat_stage
eq_flat_stage <- function(dir1, dir2, dispersal_gradient) {
  beta_res <- list()
  net_metrics_res <- list()
  
  for (i in 1:30) {
    load(paste(dir1, i, ".Rdata", sep = ""))
    stage1 <- meta_df_eq_flat %>% 
      filter(time == 1200 & N > 0) %>% 
      select(N, patch, species, time, dispersal) %>% 
      mutate(stage = "t1")
    rm(meta_df_eq_flat)
    
    load(paste(dir2, i, ".Rdata", sep = ""))
    stage2 <- meta_df_eq_flat_stage2 %>% 
      filter(time == 1200 & N > 0) %>% 
      select(N, patch, species2, time, dispersal) %>% 
      mutate(stage = "t2", patch = patch + 20) %>% 
      rename(species = species2)
    rm(meta_df_eq_flat_stage2)
    
    df <- rbind(stage1, stage2)
    beta <- list()
    net_metrics <- list()
    for (m in 1:length(dispersal_gradient)) {
      df_tep <- df %>% filter(dispersal == dispersal_gradient[m])
      webs <- bipartite::frame2webs(df_tep, varnames = c("patch", "species", "stage", "N"),
                                    emptylist = FALSE)
      beta[[m]] <- beta_temporal_metacomm(q2bnet(webs[[1]]), q2bnet(webs[[2]])) %>% 
        mutate(competition = "equal", 
               niche = "flat", 
               dispersal = dispersal_gradient[m])
      
      gamma <- sapply(webs, function(x) length(colnames(empty(x))))
      alpha <- sapply(webs, function(x) mean(rowSums(empty(x) > 0)))
      links <- sapply(webs, function(x) sum(empty(x) > 0))
      connectance <- sapply(webs, function(x) bipartite::networklevel(empty(x), index = "connectance"))
      nodf <- sapply(webs, function(x) bipartite::networklevel(empty(x), index = "NODF"))
      modularity <- sapply(webs, function(x) ifelse(length(colnames(empty(x))) > 1, 
                                                    bipartite::computeModules(empty(x))@likelihood, NA))
      
      nets_str <- data.frame(links = links, 
                             connectance = connectance,
                             nodf = nodf, 
                             modularity = modularity, 
                             gamma = gamma,
                             alpha = alpha) 
      #print(nets_str)
      nets_str$stage <- rownames(nets_str)
      res2 <- reshape2::melt(nets_str, id = "stage")
      names(res2)[2:3] <- c("structure", "value")
      net_metrics[[m]] <- res2 %>% 
        mutate(competition = "equal", 
               niche = "flat", 
               dispersal = dispersal_gradient[m])
      
    }
    
    beta_res[[i]] <- do.call("rbind", beta) %>% mutate(rep = i)
    net_metrics_res[[i]] <- do.call("rbind", net_metrics) %>% mutate(rep = i)
    print(i)
  }
  
  return(list(beta = do.call("rbind", beta_res), 
              net_metrics = do.call("rbind", net_metrics_res)))
}


##########function stab_narrow_stage
stab_narrow_stage <- function(dir1, dir2, dispersal_gradient) {
  beta_res <- list()
  net_metrics_res <- list()
  
  for (i in 1:30) {
    load(paste(dir1, i, ".Rdata", sep = ""))
    stage1 <- meta_df_stab_narrow %>% 
      filter(time == 1200 & N > 0) %>% 
      select(N, patch, species, time, dispersal) %>% 
      mutate(stage = "t1")
    rm(meta_df_stab_narrow)
    
    load(paste(dir2, i, ".Rdata", sep = ""))
    stage2 <- meta_df_stab_narrow_stage2 %>% 
      filter(time == 1200 & N > 0) %>% 
      select(N, patch, species2, time, dispersal) %>% 
      mutate(stage = "t2", patch = patch + 20) %>% 
      rename(species = species2)
    rm(meta_df_stab_narrow_stage2)
    
    df <- rbind(stage1, stage2)
    beta <- list()
    net_metrics <- list()
    for (m in 1:length(dispersal_gradient)) {
      df_tep <- df %>% filter(dispersal == dispersal_gradient[m])
      webs <- bipartite::frame2webs(df_tep, varnames = c("patch", "species", "stage", "N"),
                                    emptylist = FALSE)
      beta[[m]] <- beta_temporal_metacomm(q2bnet(webs[[1]]), q2bnet(webs[[2]])) %>% 
        mutate(competition = "stab", 
               niche = "narrow", 
               dispersal = dispersal_gradient[m])
      
      gamma <- sapply(webs, function(x) length(colnames(empty(x))))
      alpha <- sapply(webs, function(x) mean(rowSums(empty(x) > 0)))
      links <- sapply(webs, function(x) sum(empty(x) > 0))
      connectance <- sapply(webs, function(x) bipartite::networklevel(empty(x), index = "connectance"))
      nodf <- sapply(webs, function(x) bipartite::networklevel(empty(x), index = "NODF"))
      modularity <- sapply(webs, function(x) ifelse(length(colnames(empty(x))) > 1, 
                                                    bipartite::computeModules(empty(x))@likelihood, NA))
      
      nets_str <- data.frame(links = links, 
                             connectance = connectance,
                             nodf = nodf, 
                             modularity = modularity, 
                             gamma = gamma,
                             alpha = alpha) 
      #print(nets_str)
      nets_str$stage <- rownames(nets_str)
      res2 <- reshape2::melt(nets_str, id = "stage")
      names(res2)[2:3] <- c("structure", "value")
      net_metrics[[m]] <- res2 %>% 
        mutate(competition = "stab", 
               niche = "narrow", 
               dispersal = dispersal_gradient[m])
      
    }
    
    beta_res[[i]] <- do.call("rbind", beta) %>% mutate(rep = i)
    net_metrics_res[[i]] <- do.call("rbind", net_metrics) %>% mutate(rep = i)
  }
  
  return(list(beta = do.call("rbind", beta_res), 
              net_metrics = do.call("rbind", net_metrics_res)))
}


#### function stab_flat_stage
stab_flat_stage <- function(dir1, dir2, dispersal_gradient) {
  beta_res <- list()
  net_metrics_res <- list()
  
  for (i in 1:30) {
    load(paste(dir1, i, ".Rdata", sep = ""))
    stage1 <- meta_df_stab_flat %>% 
      filter(time == 1200 & N > 0) %>% 
      select(N, patch, species, time, dispersal) %>% 
      mutate(stage = "t1")
    rm(meta_df_stab_flat)
    
    load(paste(dir2, i, ".Rdata", sep = ""))
    stage2 <- meta_df_stab_flat_stage2 %>% 
      filter(time == 1200 & N > 0) %>% 
      select(N, patch, species2, time, dispersal) %>% 
      mutate(stage = "t2", patch = patch + 20) %>% 
      rename(species = species2)
    rm(meta_df_stab_flat_stage2)
    
    df <- rbind(stage1, stage2)
    beta <- list()
    net_metrics <- list()
    for (m in 1:length(dispersal_gradient)) {
      df_tep <- df %>% filter(dispersal == dispersal_gradient[m])
      webs <- bipartite::frame2webs(df_tep, varnames = c("patch", "species", "stage", "N"),
                                    emptylist = FALSE)
      beta[[m]] <- beta_temporal_metacomm(q2bnet(webs[[1]]), q2bnet(webs[[2]])) %>% 
        mutate(competition = "stab", 
               niche = "flat", 
               dispersal = dispersal_gradient[m])
      
      gamma <- sapply(webs, function(x) length(colnames(empty(x))))
      alpha <- sapply(webs, function(x) mean(rowSums(empty(x) > 0)))
      links <- sapply(webs, function(x) sum(empty(x) > 0))
      connectance <- sapply(webs, function(x) bipartite::networklevel(empty(x), index = "connectance"))
      nodf <- sapply(webs, function(x) bipartite::networklevel(empty(x), index = "NODF"))
      modularity <- sapply(webs, function(x) ifelse(length(colnames(empty(x))) > 1, 
                                                    bipartite::computeModules(empty(x))@likelihood, NA))
      
      nets_str <- data.frame(links = links, 
                             connectance = connectance,
                             nodf = nodf, 
                             modularity = modularity, 
                             gamma = gamma,
                             alpha = alpha) 
      #print(nets_str)
      nets_str$stage <- rownames(nets_str)
      res2 <- reshape2::melt(nets_str, id = "stage")
      names(res2)[2:3] <- c("structure", "value")
      net_metrics[[m]] <- res2 %>% 
        mutate(competition = "stab", 
               niche = "flat", 
               dispersal = dispersal_gradient[m])
      
    }
    
    beta_res[[i]] <- do.call("rbind", beta) %>% mutate(rep = i)
    net_metrics_res[[i]] <- do.call("rbind", net_metrics) %>% mutate(rep = i)
    print(i)
  }
  
  return(list(beta = do.call("rbind", beta_res), 
              net_metrics = do.call("rbind", net_metrics_res)))
}




#####################################################################################
## 1. calculations-----
## 1.1 equal competition scenarios-------------
t1 <- Sys.time()
t1
eq_narrow_stage <- eq_narrow_stage(dir1 = "Data/data2/meta_df_eq_narrow_", 
                                   dir2 = "Data/data2/stage2/meta_df_eq_narrow_stage2_", 
                                   dispersal_gradient = dispersal_gradient)

t2 <- Sys.time()
t2-t1

##flat niche
t1 <- Sys.time()
t1
eq_flat_stage <- eq_flat_stage(dir1 = "Data/data2/meta_df_eq_flat_", 
                                   dir2 = "Data/data2/stage2/meta_df_eq_flat_stage2_", 
                                   dispersal_gradient = dispersal_gradient)

t2 <- Sys.time()
t2-t1


## 1.2 stabilizing competition scenarios--------
##### narrow niche
t1 <- Sys.time()
t1
stab_narrow_stage <- stab_narrow_stage(dir1 = "Data/data2/meta_df_stab_narrow_", 
                                   dir2 = "Data/data2/stage2/meta_df_stab_narrow_stage2_", 
                                   dispersal_gradient = dispersal_gradient)

t2 <- Sys.time()
t2-t1

##flat niche
t1 <- Sys.time()
t1
stab_flat_stage <- stab_flat_stage(dir1 = "Data/data2/meta_df_stab_flat_", 
                                       dir2 = "Data/data2/stage2/meta_df_stab_flat_stage2_", 
                                       dispersal_gradient = dispersal_gradient)

t2 <- Sys.time()
t2-t1




## 2.make figures----
###network structure---------
fig_eq_narrow_stage_str <- ggplot(eq_narrow_stage$net_metrics %>% group_by(dispersal, structure, stage) %>% 
                                    summarise(value_m = median(na.omit(value)),
                                              value_l = quantile(na.omit(value), 0.25),
                                              value_h = quantile(na.omit(value), 0.75)), 
                                  aes(x = dispersal, y = value_m, 
                                      color = factor(stage, levels = c("t1", "t2"),
                                                     labels = c("T1", "T2")),
                                      fill = factor(stage, levels = c("t1", "t2"),
                                                    labels = c("T1", "T2")))) +
  geom_ribbon(aes(ymin = value_l, ymax = value_h), color = NA, alpha = 0.25) +
  geom_line(linewidth = 1) +
  facet_wrap(~ factor(structure, levels = c("alpha", "gamma", "links", "connectance",
                                            "nodf", "modularity"),
                      labels = c("alpha-richness", "gamma-richness", "Number of links", 
                                 "Connectance, C", "Nestedness, N", "Modularity, M")), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  scale_color_manual(values = mycol[c(3, 6)]) +
  scale_fill_manual(values = mycol[c(3, 6)]) +
  theme_bw() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Network structure") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())


fig_stab_narrow_stage_str <- ggplot(stab_narrow_stage$net_metrics %>% group_by(dispersal, structure, stage) %>% 
                                    summarise(value_m = median(na.omit(value)),
                                              value_l = quantile(na.omit(value), 0.25),
                                              value_h = quantile(na.omit(value), 0.75)), 
                                  aes(x = dispersal, y = value_m, 
                                      color = factor(stage, levels = c("t1", "t2"),
                                                     labels = c("T1", "T2")),
                                      fill = factor(stage, levels = c("t1", "t2"),
                                      labels = c("T1", "T2")))) +
  geom_ribbon(aes(ymin = value_l, ymax = value_h), color = NA, alpha = 0.25) +
  geom_line(linewidth = 1) +
  facet_wrap(~ factor(structure, levels = c("alpha", "gamma", "links", "connectance",
                                            "nodf", "modularity"),
                      labels = c("alpha-richness", "gamma-richness", "Number of links", 
                                 "Connectance, C", "Nestedness, N", "Modularity, M")), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  scale_color_manual(values = mycol[c(3, 6)]) +
  scale_fill_manual(values = mycol[c(3, 6)]) +
  theme_bw() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Network structure") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())


fig_eq_narrow_stage_str + fig_stab_narrow_stage_str + plot_layout(ncol = 1) +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(', tag_sep = '', 
                  tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 12,family = "serif", hjust = 0, vjust = 0))

ggsave("./Figs/R1/Figure_stage_netstr.png", width = 6*1.2, height = 16*0.55)
ggsave("./Figs/R1/Figure_stage_netstr.pdf", width = 9*1.2, height = 16*0.55)

### beta diversity----------
eq_stab_stage_narrow_beta <- rbind(eq_narrow_stage$beta %>% group_by(dispersal, Component) %>% 
                                     summarise(beta_m = median(beta),
                                               beta_l = quantile(beta, 0.25),
                                               beta_h = quantile(beta, 0.75)) %>% 
                                     mutate(niche = "narrow", competition = "equal competition"),
                                   stab_narrow_stage$beta %>% group_by(dispersal, Component) %>% 
                                     summarise(beta_m = median(beta),
                                               beta_l = quantile(beta, 0.25),
                                               beta_h = quantile(beta, 0.75)) %>% 
                                     mutate(niche = "narrow", competition = "stabilizing competition"))


ggplot(eq_stab_stage_narrow_beta, 
       aes(x = dispersal, y = beta_m, color = factor(competition),
           fill = factor(competition))) +
  geom_ribbon(aes(ymin = beta_l, ymax = beta_h), 
              color = NA, alpha = 0.45) +
  geom_line(size = 1) +
  facet_wrap(~ factor(Component, levels = c("beta_temporal", "beta_colonization", "beta_extinction",
                                            "beta_Local", "beta_Local_colonization", "beta_Local_extinction",
                                            "beta_Regional", "beta_Regional_colonization", "beta_Regional_extinction",
                                            "beta_Landscape", "beta_Landscape_gain", "beta_Landscape_loss", 
                                            "beta_RL", "beta_RL_colonization","beta_RL_extinction"),
                      labels = c("Total beta", "Total colonization", "Total extinction",
                                 "Local", "Local colonization", "Local extinction",
                                 "Regional", "Regional colonization", " Regional extinction",
                                 "Landscape", "Patch gain ", "Patch loss", 
                                 "Reginal landscape (RL)", "RL colonization",
                                 "RL extinction")), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  scale_color_manual(values = mycol[c(3, 6)]) +
  scale_fill_manual(values = mycol[c(3, 6)]) +
  theme_bw() +
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Temporal beta diversity of species-patch links") +
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())

ggsave("./Figs/R1/eq_stab_narrow_stage_beta.png",
       width = 8.2, height = 8.2, units = "in", dpi = 300)

ggsave("./Figs/R1/eq_stab_narrow_stage_beta.pdf",
       width = 8.2, height = 8.2, units = "in")

################################################################################
fig_eq_flat_stage_str <- ggplot(eq_flat_stage$net_metrics %>% group_by(dispersal, structure, stage) %>% 
                                    summarise(value_m = median(na.omit(value)),
                                              value_l = quantile(na.omit(value), 0.25),
                                              value_h = quantile(na.omit(value), 0.75)), 
                                  aes(x = dispersal, y = value_m, 
                                      color = factor(stage, levels = c("t1", "t2"),
                                                     labels = c("T1", "T2")),
                                      fill = factor(stage, levels = c("t1", "t2"),
                                                    labels = c("T1", "T2")))) +
  geom_ribbon(aes(ymin = value_l, ymax = value_h), color = NA, alpha = 0.25) +
  geom_line(linewidth = 1) +
  facet_wrap(~ factor(structure, levels = c("alpha", "gamma", "links", "connectance",
                                            "nodf", "modularity"),
                      labels = c("alpha-richness", "gamma-richness", "Number of links", 
                                 "Connectance, C", "Nestedness, N", "Modularity, M")), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  scale_color_manual(values = mycol[c(3, 6)]) +
  scale_fill_manual(values = mycol[c(3, 6)]) +
  theme_bw() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Network structure") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())


fig_stab_flat_stage_str <- ggplot(stab_flat_stage$net_metrics %>% group_by(dispersal, structure, stage) %>% 
                                      summarise(value_m = median(na.omit(value)),
                                                value_l = quantile(na.omit(value), 0.25),
                                                value_h = quantile(na.omit(value), 0.75)), 
                                    aes(x = dispersal, y = value_m, 
                                        color = factor(stage, levels = c("t1", "t2"),
                                                       labels = c("T1", "T2")),
                                        fill = factor(stage, levels = c("t1", "t2"),
                                                      labels = c("T1", "T2")))) +
  geom_ribbon(aes(ymin = value_l, ymax = value_h), color = NA, alpha = 0.25) +
  geom_line(linewidth = 1) +
  facet_wrap(~ factor(structure, levels = c("alpha", "gamma", "links", "connectance",
                                            "nodf", "modularity"),
                      labels = c("alpha-richness", "gamma-richness", "Number of links", 
                                 "Connectance, C", "Nestedness, N", "Modularity, M")), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  scale_color_manual(values = mycol[c(3, 6)]) +
  scale_fill_manual(values = mycol[c(3, 6)]) +
  theme_bw() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Network structure") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())


fig_eq_narrow_stage_str + fig_eq_flat_stage_str + plot_layout(ncol = 1) +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(', tag_sep = '', 
                  tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 12,family = "serif", hjust = 0, vjust = 0))

ggsave("./Figs/R1/Figure_eq_stage_netstr.png", width = 6*1.2, height = 16*0.55)
ggsave("./Figs/R1/Figure_eq_stage_netstr.pdf", width = 9*1.2, height = 16*0.55)


fig_stab_narrow_stage_str + fig_stab_flat_stage_str + plot_layout(ncol = 1) +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(', tag_sep = '', 
                  tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 12,family = "serif", hjust = 0, vjust = 0))

ggsave("./Figs/R1/Figure_stab_stage_netstr.png", width = 6*1.2, height = 16*0.55)
ggsave("./Figs/R1/Figure_stab_stage_netstr.pdf", width = 9*1.2, height = 16*0.55)

##################################################################################

eq_stage_beta <- rbind(eq_narrow_stage$beta %>% group_by(dispersal, Component) %>% 
                                     summarise(beta_m = median(beta),
                                               beta_l = quantile(beta, 0.25),
                                               beta_h = quantile(beta, 0.75)) %>% 
                                     mutate(niche = "Narrow niche", competition = "equal competition"),
                                   eq_flat_stage$beta %>% group_by(dispersal, Component) %>% 
                                     summarise(beta_m = median(beta),
                                               beta_l = quantile(beta, 0.25),
                                               beta_h = quantile(beta, 0.75)) %>% 
                                     mutate(niche = "Wide niche", competition = "equal competition"))


ggplot(eq_stage_beta, 
       aes(x = dispersal, y = beta_m, color = factor(niche),
           fill = factor(niche))) +
  geom_ribbon(aes(ymin = beta_l, ymax = beta_h), 
              color = NA, alpha = 0.45) +
  geom_line(size = 1) +
  facet_wrap(~ factor(Component, levels = c("beta_temporal", "beta_colonization", "beta_extinction",
                                            "beta_Local", "beta_Local_colonization", "beta_Local_extinction",
                                            "beta_Regional", "beta_Regional_colonization", "beta_Regional_extinction",
                                            "beta_Landscape", "beta_Landscape_gain", "beta_Landscape_loss", 
                                            "beta_RL", "beta_RL_colonization","beta_RL_extinction"),
                      labels = c("Total beta", "Total colonization", "Total extinction",
                                 "Local", "Local colonization", "Local extinction",
                                 "Regional", "Regional colonization", " Regional extinction",
                                 "Landscape", "Patch gain ", "Patch loss", 
                                 "Reginal landscape (RL)", "RL colonization",
                                 "RL extinction")), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  scale_color_manual(values = mycol[c(3, 6)]) +
  scale_fill_manual(values = mycol[c(3, 6)]) +
  theme_bw() +
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Temporal beta diversity of species-patch links") +
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())

ggsave("./Figs/R1/eq_stage_beta.png",
       width = 8.2, height = 8.2, units = "in", dpi = 300)

ggsave("./Figs/R1/eq_stage_beta.pdf",
       width = 8.2, height = 8.2, units = "in")

################
stab_stage_beta <- rbind(stab_narrow_stage$beta %>% group_by(dispersal, Component) %>% 
                         summarise(beta_m = median(beta),
                                   beta_l = quantile(beta, 0.25),
                                   beta_h = quantile(beta, 0.75)) %>% 
                         mutate(niche = "Narrow niche", competition = "stabilizing competition"),
                       stab_flat_stage$beta %>% group_by(dispersal, Component) %>% 
                         summarise(beta_m = median(beta),
                                   beta_l = quantile(beta, 0.25),
                                   beta_h = quantile(beta, 0.75)) %>% 
                         mutate(niche = "Wide niche", competition = "stabilizing competition"))


ggplot(stab_stage_beta, 
       aes(x = dispersal, y = beta_m, color = factor(niche),
           fill = factor(niche))) +
  geom_ribbon(aes(ymin = beta_l, ymax = beta_h), 
              color = NA, alpha = 0.45) +
  geom_line(size = 1) +
  facet_wrap(~ factor(Component, levels = c("beta_temporal", "beta_colonization", "beta_extinction",
                                            "beta_Local", "beta_Local_colonization", "beta_Local_extinction",
                                            "beta_Regional", "beta_Regional_colonization", "beta_Regional_extinction",
                                            "beta_Landscape", "beta_Landscape_gain", "beta_Landscape_loss", 
                                            "beta_RL", "beta_RL_colonization","beta_RL_extinction"),
                      labels = c("Total beta", "Total colonization", "Total extinction",
                                 "Local", "Local colonization", "Local extinction",
                                 "Regional", "Regional colonization", " Regional extinction",
                                 "Landscape", "Patch gain ", "Patch loss", 
                                 "Reginal landscape (RL)", "RL colonization",
                                 "RL extinction")), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  scale_color_manual(values = mycol[c(3, 6)]) +
  scale_fill_manual(values = mycol[c(3, 6)]) +
  theme_bw() +
  xlab(expression(paste("Dispersal (",italic(a[i]),")",sep = ""))) + 
  ylab("Temporal beta diversity of species-patch links") +
  theme_bw()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())

ggsave("./Figs/R1/stab_stage_beta.png",
       width = 8.2, height = 8.2, units = "in", dpi = 300)

ggsave("./Figs/R1/stab_stage_beta.pdf",
       width = 8.2, height = 8.2, units = "in")

