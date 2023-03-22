###################################################################################
####### sensitivity of patches for Horvath data####################################
###################################################################################
net_sampler <- function(net, pro = seq(0.1,1, by = 0.1), nsim = 10){
  n_patch <- nrow(net)
  
  net_metrics <- list()
  
  for (i in 1:length(pro)) {
    n_sample <- round(pro[i]*n_patch)
    
    webs <- list()
    for (sim in 1:nsim) {
      web_tep <- bipartite::empty(net[sample(1:n_patch, n_sample),])
      names(web_tep) <- sim
      webs[[sim]] <- web_tep
    }
    
    gamma <- sapply(webs, function(x) length(colnames(x)))
    alpha <- sapply(webs, function(x) mean(rowSums(x > 0)))
    links <- sapply(webs, function(x) sum(x > 0))
    connectance <- sapply(webs, function(x) bipartite::networklevel(x, index = "connectance"))
    nodf <- sapply(webs, function(x) bipartite::networklevel(x, index = "NODF"))
    modularity <- sapply(webs, function(x) ifelse(length(colnames(x)) > 1, 
                                                  bipartite::computeModules(x)@likelihood, NA))
    
    nets_str <- data.frame(links = links, 
                           connectance = connectance,
                           nodf = nodf, 
                           modularity = modularity, 
                           gamma = gamma,
                           alpha = alpha) 
    #print(nets_str)
    nets_str$sim <- rownames(nets_str)
    res2 <- reshape2::melt(nets_str, id = "sim")
    names(res2)[2:3] <- c("structure", "value")
    net_metrics[[i]] <- res2 %>% mutate(effort = pro[i])
  }
  
  net_metrics_res <- do.call("rbind", net_metrics)
  net_metrics_res
}



beta_sampler <- function(net1, net2, pro = seq(0.1, 1, by = 0.1), nsim = 10) {
  # all dimension of net1 and net2 should be same!!!
  n_patch <- nrow(net1)
  beta <- list()
  for (i in 1:length(pro)) {
    
    n_sample <- round(pro[i]*n_patch)
    beta_sim <- list()
    
    for (sim in 1:nsim){
      web_tep1 <- net1[sample(1:n_patch, n_sample),]
      web_tep2 <- net2[rownames(web_tep1),]
      beta_sim[[sim]] <- beta_temporal_metacomm(web_tep1, web_tep2) %>% 
        mutate(sim = sim) 
    }
    beta[[i]] <- do.call("rbind", beta_sim) %>% mutate(effort = pro[i])
    print(i/length(pro))
  }
  beta_res <- do.call("rbind", beta)
  beta_res
}



## data import and check
Horvath <- read.xlsx("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Data/Horvath_et_al_2019_ Ecology_Letters_community_env_data.xlsx", 
                     sheet = 2)
names(Horvath)
unique(Horvath$Year)

Horvath_1957 <- Horvath %>% 
  filter(Year == 1957) %>% 
  dplyr::select(c(2,9:77))

Horvath_1957 <- right.matrix(Horvath_1957)
visweb(Horvath_1957)
dim(Horvath_1957)

Horvath_2010 <- Horvath %>% 
  filter(Year == 2010) %>% 
  dplyr::select(c(2,9:77))

Horvath_2010 <- right.matrix(Horvath_2010)
visweb(Horvath_2010)

Horvath_2010 <- Horvath_2010[rownames(Horvath_1957),] # make rownames equal

all(rownames(Horvath_1957)==rownames(Horvath_2010)) # should be true


all(colSums(Horvath_1957)+colSums(Horvath_2010)>0)


dim(Horvath_1957) # should be 53 patches


Horvath_nets <- list(Horvath_1957 = Horvath_1957,
                     Horvath_2020 = Horvath_2010)

### simulation
### network structure
net_str_horvath_1957 <- net_sampler(net = Horvath_1957, nsim = 100)

###temporal beta diversity
net_beta_horvath <- beta_sampler(net1 = Horvath_1957, net2 = Horvath_2010, nsim = 100)

unique(net_beta_horvath$Component)
### make figures

ggplot(net_str_horvath_1957, aes(x = as.numeric(effort), y = value)) +
  geom_point(size = 2, alpha = 0.15, color = mycol[3]) +
  geom_smooth(method = "loess", se = F,size = 1.2, color = mycol[6]) +
  facet_wrap(~ factor(structure, levels = c("alpha","gamma", "links", 
                                            "connectance", "nodf", "modularity"),
                      labels = c("alpha-richness","gamma-richness", "Number of links", 
                                 "Connectance, C", "Nestedness, NODF", "Modularity, M")), 
             ncol = 3, scales = "free") + 
  theme_bw() +
  xlab("Sampling effort of patches") + 
  ylab("Species richness and network structure") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 12,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())

ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/horvath_str_sensitivity.pdf",
       width = 7.8, height = 5.8, units = "in")

ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/horvath_str_sensitivity.png",
       width = 7.8, height = 5.8, units = "in", dpi = 600)

### beta vs. patch number
ggplot(net_beta_horvath %>% filter(Component == "beta_temporal" | 
                                     Component == "beta_extinction" |
                                     Component == "beta_colonization" |
                                     Component == "beta_Local" |
                                     Component == "beta_Regional" |
                                     Component == "beta_Landscape" |
                                     Component == "beta_RL" ), 
       aes(x = as.numeric(effort), y = beta)) +
  geom_point(size = 2, alpha = 0.1, color = mycol[3]) +
  geom_smooth(method = "loess", se = F,size = 1.2, color = mycol[6]) +
  facet_wrap(~ factor(Component, levels = c("beta_temporal", "beta_colonization", "beta_extinction", "",
                                            "beta_Local", "beta_Regional", 
                                            "beta_Landscape", "beta_RL")), 
             ncol = 4, drop = FALSE, 
             scales = "free") + 
  #scale_color_manual(values = mycol[1:3]) +
  theme_bw() +
  xlab("Sampling effort of patches") + 
  ylab("Beta diversity") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 12,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())
ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/horvath_beta_sensitivity_1.pdf",
       width = 10, height = 5.8, units = "in")


ggplot(net_beta_horvath %>% filter(Component != "beta_Landscape_gain" &
                                     Component != "beta_RL_colonization"), 
       aes(x = as.numeric(effort), y = beta)) +
  geom_point(size = 2, alpha = 0.1, color = mycol[3]) +
  geom_smooth(method = "loess", se = F,size = 1.2, color = mycol[6]) +
  facet_wrap(~ factor(Component, levels = c("beta_temporal", "beta_colonization", "beta_extinction",
                                            "beta_Local", "beta_Local_colonization", "beta_Local_extinction",
                                            "beta_Regional", "beta_Regional_colonization", "beta_Regional_extinction",
                                            "beta_Landscape", " ", "beta_Landscape_loss", 
                                            "beta_RL", "","beta_RL_extinction"),
                      labels = c("Total beta", "Total colonization", "Total extinction",
                                 "Local", "Local colonization", "Local extinction",
                                 "Regional", "Regional colonization", " Regional extinction",
                                 "Landscape", "Patch gain ", "Patch loss", 
                                 "Reginal landscape", "Reginal landscape colonization",
                                 "Reginal landscape extinction")), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  #scale_color_manual(values = mycol[1:3]) +
  theme_bw() +
  xlab("Sampling effort of patches") + 
  ylab("Beta diversity") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())

ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/horvath_beta_sensitivity_2.pdf",
       width = 7.2, height = 8.2, units = "in")
ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/horvath_beta_sensitivity_2.png",
       width = 7.2, height = 8.2, units = "in", dpi = 300)



ggplot(net_beta_horvath, 
       aes(x = as.numeric(effort), y = beta)) +
  geom_point(size = 2, alpha = 0.1, color = mycol[3]) +
  geom_smooth(method = "loess", se = F,size = 1.2, color = mycol[6]) +
  facet_wrap(~ factor(Component), 
             ncol = 3, drop = FALSE,
             scales = "free") + 
  #scale_color_manual(values = mycol[1:3]) +
  theme_bw() +
  xlab("Sampling effort of patches") + 
  ylab("Beta diversity") +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.y = element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12,family = "serif"),
        strip.text = element_text(size = 11,family = "serif"),
        legend.text = element_text(size = 12,family = "serif"),
        legend.title = element_blank())

