## Rcode of the case studies
library(openxlsx)
library(tydiverse)
library(bipartite)
library(igraph)
library(lme4)
library(MuMIn)

# case study 1 Hoja mammal data-------------
Hoja_11 <- read.xlsx("C:/Dataset/Temporal_metacommunity/temporal_beta_metacommunity/Data/Hoja.xlsx",
                     sheet = 2,rowNames = T)

Hoja_17 <- read.xlsx("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Data/Hoja.xlsx",
                     sheet = 3,rowNames = T)

networklevel(Hoja_11,index = "connectance")
## connectance 
## 0.2928571

networklevel(Hoja_17,index = "connectance")
## connectance 
## 0.5571429

networklevel(Hoja_11,index = "NODF") 
# NODF 
# 49.64975

networklevel(Hoja_17,index = "NODF")
# NODF 
# 76.18112

# here we can see that the nestedness is increasing

nodf_Hoja_obs <- sapply(list(Hoja_11, Hoja_17), 
                        networklevel, index = "NODF")
nodf_Hoja_obs

null_nets_11 <- simulate(vegan::nullmodel(Hoja_11, method = "quasiswap"), 
                         nsim = 1000)

null_nets_17 <- simulate(vegan::nullmodel(Hoja_17, method = "quasiswap"), 
                         nsim = 1000)

nodf_nulls_11 <- apply(null_nets_11, 3, networklevel, index= "NODF")

hist(nodf_nulls_11)

z_11 <- (nodf_Hoja_obs[1] - mean(nodf_nulls_11)) / sd(nodf_nulls_11)
z_11

# NODF 
# -0.059814


nodf_nulls_17 <- apply(null_nets_17, 3, networklevel, index= "NODF")

hist(nodf_nulls_17)

z_17 <- (nodf_Hoja_obs[2] - mean(nodf_nulls_17)) / sd(nodf_nulls_17)
z_17
# NODF 
# 1.274885

# species contribute to nestedness---
cnodf_11 <- nestedcontribution(Hoja_11, nsimul = 999)

cnodf_17 <- nestedcontribution(Hoja_17, nsimul = 999)

# delta degree and delta condf
degree_sp <- merge(specieslevel(Hoja_11,index = "degree")$`higher level` %>% 
                     mutate(species = rownames(specieslevel(Hoja_11, index = "degree")$`higher level`),
                            year = "2011") %>% 
                     rename(degree_2011 = degree),
                   specieslevel(Hoja_17, index = "degree")$`higher level` %>% 
                     mutate(species = rownames(specieslevel(Hoja_17, index = "degree")$`higher level`),
                            year = "2017")%>% 
                     rename(degree_2017 = degree),
                   by = "species") %>% 
  mutate(delta_degree = degree_2017 - degree_2011) %>% 
  ggplot(.,aes(x = reorder(species, -delta_degree), y = delta_degree)) +
  geom_bar(stat = "identity", fill = "grey80") +
  geom_hline(yintercept = 0) +
  ylab("Delta degree") + xlab("Species") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.title.y = element_text(size = 11,family = "serif"),
        axis.text.y = element_text(size = 11, family = "serif"),
        axis.text.x = element_text(size = 11, family = "serif",
                                   angle = 90,vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 11,family = "serif"))

degree_sp


cnodf_sp <- merge(cnodf_11$`higher level` %>% 
                    mutate(species = rownames(cnodf_11$`higher level`),
                           year = "2011") %>% 
                    rename(cnodf_2011 = nestedcontribution),
                  cnodf_17$`higher level` %>% 
                    mutate(species = rownames(cnodf_17$`higher level`),
                           year = "2017")%>% 
                    rename(cnodf_2017 = nestedcontribution),
                  by = "species") %>% 
  mutate(delta_cnodf = cnodf_2017 - cnodf_2011) %>% 
  ggplot(.,aes(x = reorder(species, -delta_cnodf), y = delta_cnodf)) +
  geom_bar(stat = "identity", fill = "grey80") +
  geom_hline(yintercept = 0)+
  ylab("Delta cnodf") + xlab("Species") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.title.y = element_text(size = 11,family = "serif"),
        axis.text.y = element_text(size = 11, family = "serif"),
        axis.text.x = element_text(size = 11, family = "serif",
                                   angle = 90,vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 11,family = "serif"))
cnodf_sp

degree_patch <- merge(specieslevel(Hoja_11,index = "degree")$`lower level` %>% 
                        mutate(patch = rownames(specieslevel(Hoja_11,index = "degree")$`lower level`),
                               year = "2011") %>% 
                        rename(degree_2011 = degree),
                      specieslevel(Hoja_17,index = "degree")$`lower level` %>% 
                        mutate(patch = rownames(specieslevel(Hoja_17,index = "degree")$`lower level`),
                               year = "2017")%>% 
                        rename(degree_2017 = degree),
                      by = "patch") %>% 
  mutate(delta_degree = degree_2017 - degree_2011) %>% 
  ggplot(.,aes(x = reorder(patch, -delta_degree), y = delta_degree)) +
  geom_bar(stat = "identity", fill = "grey80") +
  geom_hline(yintercept = 0) +
  ylab("Delta degree") + xlab("Patch") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.title.y=element_text(size = 11,family = "serif"),
        axis.text.y = element_text(size = 11, family = "serif"),
        axis.text.x = element_text(size = 11, family = "serif",
                                   angle = 90,vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 11,family = "serif"))
degree_patch




cnodf_patch <- merge(cnodf_11$`lower level` %>% 
                       mutate(patch = rownames(cnodf_11$`lower level`),
                              year = "2011") %>% 
                       rename(cnodf_2011 = nestedcontribution),
                     cnodf_17$`lower level` %>% 
                       mutate(patch = rownames(cnodf_17$`lower level`),
                              year = "2017")%>% 
                       rename(cnodf_2017 = nestedcontribution),
                     by = "patch") %>% 
  mutate(delta_cnodf = cnodf_2017-cnodf_2011) %>% 
  ggplot(.,aes(x = reorder(patch, -delta_cnodf), y = delta_cnodf)) +
  geom_bar(stat = "identity", fill = "grey80")+
  geom_hline(yintercept = 0) +
  ylab("Delta cnodf") + xlab("Patch") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.title.y = element_text(size = 11,family = "serif"),
        axis.text.y = element_text(size = 11, family = "serif"),
        axis.text.x = element_text(size = 11, family = "serif",
                                   angle = 90,vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 11,family = "serif"))
cnodf_patch


degree_sp + cnodf_sp + degree_patch + cnodf_patch + plot_spacer() +
  Hoja_beta_plot + plot_layout(ncol = 2)


ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/Hoja_beta_netstr_plot.pdf",
       height = 297*0.6, width = 210, units = "mm")





delta_species <- merge(specieslevel(Hoja_11,index = "degree")$`higher level` %>% 
                         mutate(species = rownames(specieslevel(Hoja_11,index = "degree")$`higher level`),
                                year = "2011") %>% 
                         rename(degree_2011 = degree),
                       specieslevel(Hoja_17,index = "degree")$`higher level` %>% 
                         mutate(species = rownames(specieslevel(Hoja_17,index = "degree")$`higher level`),
                                year = "2017")%>% 
                         rename(degree_2017 = degree),
                       by = "species") %>% 
  mutate(delta_degree = degree_2017-degree_2011) %>% 
  left_join(.,merge(cnodf_11$`higher level` %>% 
                      mutate(species = rownames(cnodf_11$`higher level`),
                             year = "2011") %>% 
                      rename(cnodf_2011 = nestedcontribution),
                    cnodf_17$`higher level` %>% 
                      mutate(species = rownames(cnodf_17$`higher level`),
                             year = "2017")%>% 
                      rename(cnodf_2017 = nestedcontribution),
                    by = "species") %>% 
              mutate(delta_cnodf = cnodf_2017-cnodf_2011), id = "patch")


ggplot(delta_species, aes(x = delta_degree, y = delta_cnodf)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess")

summary(mgcv::gam(delta_cnodf ~ s(delta_degree), data = delta_species))


delta_patch <- read.xlsx("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Data/Hoja.xlsx", 
                         sheet = 5) %>% 
  left_join(.,merge(cnodf_11$`lower level` %>% 
                      mutate(patch = rownames(cnodf_11$`lower level`),
                             year = "2011") %>% 
                      rename(cnodf_2011 = nestedcontribution),
                    cnodf_17$`lower level` %>% 
                      mutate(patch = rownames(cnodf_17$`lower level`),
                             year = "2017")%>% 
                      rename(cnodf_2017 = nestedcontribution),
                    by = "patch")%>% 
              mutate(delta_cnodf = cnodf_2017-cnodf_2011), id = "patch") %>% 
  left_join(.,merge(specieslevel(Hoja_11,index = "degree")$`lower level` %>% 
                      mutate(patch = rownames(specieslevel(Hoja_11,index = "degree")$`lower level`),
                             year = "2011") %>% 
                      rename(degree_2011 = degree),
                    specieslevel(Hoja_17,index = "degree")$`lower level` %>% 
                      mutate(patch = rownames(specieslevel(Hoja_17,index = "degree")$`lower level`),
                             year = "2017")%>% 
                      rename(degree_2017 = degree),
                    by = "patch")%>% 
              mutate(delta_degree = degree_2017-degree_2011), id = "patch")


names(delta_patch)
delta_patch <- delta_patch[,-c(38,40)]

pairs(delta_patch[,-1])

ggplot(delta_patch, aes(x = delta_degree, y = delta_cnodf)) +
  geom_point(size = 2) +
  geom_smooth(method = "gam")

summary(mgcv::gam(delta_cnodf ~ s(delta_degree), data = delta_patch))







rbind(cnodf_11$`higher level` %>% 
        mutate(species = rownames(cnodf_11$`higher level`),
               year = "2011"),
      cnodf_17$`higher level` %>% 
        mutate(species = rownames(cnodf_17$`higher level`),
               year = "2017")) %>% 
  ggplot(.,aes(x = species, y = nestedcontribution, group = species, fill = year))+
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"))

rbind(specieslevel(Hoja_11,index = "degree")$`higher level` %>% 
        mutate(species = rownames(specieslevel(Hoja_11,index = "degree")$`higher level`),
               year = "2011"),
      specieslevel(Hoja_17,index = "degree")$`higher level` %>% 
        mutate(species = rownames(specieslevel(Hoja_17,index = "degree")$`higher level`),
               year = "2017")) %>% 
  ggplot(.,aes(x = species, y = degree, group = species, fill = year)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"))


rbind(cnodf_11$`lower level` %>% 
        mutate(patch = rownames(cnodf_11$`lower level`),
               year = "2011"),
      cnodf_17$`lower level` %>% 
        mutate(patch = rownames(cnodf_17$`lower level`),
               year = "2017")) %>% 
  ggplot(.,aes(x = patch, y = nestedcontribution, group = patch, fill = year))+
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"))


rbind(specieslevel(Hoja_11,index = "degree")$`lower level` %>% 
        mutate(patch = rownames(specieslevel(Hoja_11,index = "degree")$`lower level`),
               year = "2011"),
      specieslevel(Hoja_17,index = "degree")$`lower level` %>% 
        mutate(patch = rownames(specieslevel(Hoja_17,index = "degree")$`lower level`),
               year = "2017")) %>% 
  ggplot(.,aes(x = patch, y = degree, group = patch, fill = year)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"))



### modularity
M_Hoja_obs <- sapply(list(Hoja_11, Hoja_17), 
                     function(x)computeModules(x)@likelihood)
M_Hoja_obs

#[1] 0.2561306 0.1651610



plotweb(Hoja_11, labsize = 0.5, method = "normal", col.high = "#8491B4B2",
        col.low = "#00A087B2", col.interaction = "grey80")

plotweb(Hoja_17, labsize = 0.5, method = "normal", col.high = "#8491B4B2",
        col.low = "#00A087B2", col.interaction = "grey80")

Hoja_beta <- beta_temporal_metacomm(Hoja_11, Hoja_17)
Hoja_beta

##         Component      beta
##        beta_Local 0.6054688
##     beta_Regional 0.0000000
##    beta_Landscape 0.0000000
##           beta_RL 0.0000000
##     beta_temporal 0.6054688
##   beta_extinction 0.0859375
## beta_colonization 0.5195312

Hoja_beta_plot <- ggplot(Hoja_beta,aes(x = Component, y = beta, fill = Component))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("β diversity") + xlab("")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 12, family = "serif"))

Hoja_beta_plot
ggsave("C:/Dataset/Temporal_metacommunity/Figs/Hoja_beta_plot.pdf", plot = Hoja_beta_plot,
       height = 395/3, width = 421/3, units = "mm")

ggsave("C:/Dataset/Temporal_metacommunity/Figs/Hoja_beta_plot.png", plot = Hoja_beta_plot,
       height = 395/3, width = 421/3, units = "mm")

# case study 2: Zooplankton data----

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

Horvath_2010 <- Horvath_2010[rownames(Horvath_1957),]
all(rownames(Horvath_1957)==rownames(Horvath_2010))


all(colSums(Horvath_1957)+colSums(Horvath_2010)>0)

dim(Horvath_1957)

Horvath_nets <- list(Horvath_1957 = Horvath_1957,
                     Horvath_2020 = Horvath_2010)


sapply(Horvath_nets, function(x)networklevel(empty(x), 
                                             index = c("connectance", "number of species", "NODF")))

64-43

# Horvath_1957
# connectance             0.1751179
# NODF                   32.9008512
# number.of.species.HL   64.0000000
# number.of.species.LL   53.0000000

# Horvath_2020
# connectance             0.1841085
# NODF                   39.2164904
# number.of.species.HL   43.0000000
# number.of.species.LL   24.0000000

Horvath_nodf_obs <- sapply(Horvath_nets, function(x)networklevel(empty(x), 
                                                                 index = "NODF"))

# Horvath_1957.NODF Horvath_2020.NODF 
# 32.90085          39.21649

null_Horvath_1957 <- simulate(vegan::nullmodel(empty(Horvath_1957), method = "quasiswap"), 
                              nsim = 1000)

null_Horvath_2010 <- simulate(vegan::nullmodel(empty(Horvath_2010), method = "quasiswap"), 
                              nsim = 1000)

null_Horvath_1957_nodf <- apply(null_Horvath_1957, 3, networklevel, index= "NODF")

Horvath_1957_nodf_z <- (Horvath_nodf_obs[1]-mean(null_Horvath_1957_nodf))/sd(null_Horvath_1957_nodf)
Horvath_1957_nodf_z
# -0.07393968

plot(density(null_Horvath_1957_nodf))
abline(v = Horvath_nodf_obs[1], col = "red")

null_Horvath_2010_nodf <- apply(null_Horvath_2010, 3, networklevel, index= "NODF")

Horvath_2010_nodf_z <- (Horvath_nodf_obs[2]-mean(null_Horvath_2010_nodf))/sd(null_Horvath_2010_nodf)
Horvath_2010_nodf_z
#0.8231277

plot(density(null_Horvath_2010_nodf))
abline(v = Horvath_nodf_obs[2], col = "red")


mudule_1957 <- computeModules(empty(Horvath_1957))
plotModuleWeb(mudule_1957)

mudule_2010 <- computeModules(empty(Horvath_2010))
plotModuleWeb(mudule_2010)

set.seed(123)
Horvath_m_obs <- sapply(Horvath_nets, function(x)computeModules(empty(x))@likelihood)
# Horvath_1957 Horvath_2020 
# 0.2963813    0.3372022

null_Horvath_1957_m <- apply(null_Horvath_1957[,,1:100], 3, function(x)computeModules(x)@likelihood)

Horvath_1957_m_z <- (Horvath_m_obs[1]-mean(null_Horvath_1957_m))/sd(null_Horvath_1957_m)
Horvath_1957_m_z
#14.53441 

null_Horvath_2010_m <- apply(null_Horvath_2010[,,1:100], 3, function(x)computeModules(x)@likelihood)

Horvath_2010_m_z <- (Horvath_m_obs[2]-mean(null_Horvath_2010_m))/sd(null_Horvath_2010_m)
Horvath_2010_m_z


z_score_net <- function(net){
  nodf_obs <- bipartite::networklevel(net, index = "NODF")
  M_obs <- bipartite::computeModules(net)@likelihood
  nulls <- simulate(vegan::nullmodel(net, method = "quasiswap"), nsim = 1000)
  null_nodf <- apply(nulls, 3, networklevel, index= "NODF")
  z_nodf <- (nodf_obs-mean(null_nodf))/sd(null_nodf)
  
  null_M <-apply(nulls, 3, function(x)computeModules(x)@likelihood)
  z_M <- (M_obs-mean(M_nodf))/sd(null_M)
  
  res <- c(nodf_obs, z_nodf, M_obs, z_M)
  res
  
}




Horvath_beta <- beta_temporal_metacomm(Horvath_1957, Horvath_2010)
Horvath_beta

#           Component       beta
# 1        beta_Local 0.33991537
# 2     beta_Regional 0.07334274
# 3    beta_Landscape 0.35119887
# 4           beta_RL 0.12976023
# 5     beta_temporal 0.89421721
# 6   beta_extinction 0.73201693
# 7 beta_colonization 0.16220028

# 0.73201693/0.16220028 = 4.513044
# 0.73201693/0.89421721 = 0.818612

Horvath_beta %>% 
  mutate(proportion = beta/0.89421721)

#           Component       beta proportion
# 1        beta_Local 0.33991537 0.38012618
# 2     beta_Regional 0.07334274 0.08201893
# 3    beta_Landscape 0.35119887 0.39274448
# 4           beta_RL 0.12976023 0.14511041
# 5     beta_temporal 0.89421721 1.00000000
# 6   beta_extinction 0.73201693 0.81861198
# 7 beta_colonization 0.16220028 0.18138801


Horvath_beta_plot <- ggplot(Horvath_beta,aes(x = Component, y = beta, fill = Component))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("β diversity") + xlab("")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 12, family = "serif"))
Horvath_beta_plot


ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/Hovath_beta_plot.pdf",
       plot = Horvath_beta_plot,
       height = 395/3, width = 421/3, units = "mm")

op <- par(mfrow = c(2,1))
plotweb(Horvath_1957, method = "norm",empty = F)
plotweb(Horvath_2010, method = "norm",empty = F)
par(op)

plotweb(Horvath_2010, empty = T)
                    

# case study 3 Lindholm 2019 warter plant data-------------

# a function to make a dataframe to matrix
right.matrix <- function(x) {
  m<-as.matrix(x[, -1])
  rownames(m) <- x[, 1]
  m
}

Lindholm <- read.xlsx("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Data/Lindholm_2019.xlsx", sheet = 2)
unique(Lindholm$Year)

Lindholm_env <- read.csv2("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Data/Lindholm_2019_Env.csv",
                          sep = ";")
head(Lindholm_env)

write.xlsx(Lindholm_env, file = "C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Data/Lindholm_2019_Env.xlsx")

Lindholm_env <- read.xlsx("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Data/Lindholm_2019_Env.xlsx", 
                          sheet = 1, rowNames = T)


names(Lindholm_env)

Lindholm_1940 <- Lindholm %>% filter(Year == 1940) %>% 
  dplyr::select(c(1,3:68))

Lindholm_1940 <- right.matrix(Lindholm_1940)

Lindholm_1970 <- Lindholm %>% filter(Year == 1970) %>% 
  dplyr::select(c(1,3:68))

Lindholm_1970 <- right.matrix(Lindholm_1970)

Lindholm_1990 <- Lindholm %>% filter(Year == 1990) %>% 
  dplyr::select(c(1,3:68))
Lindholm_1990 <- right.matrix(Lindholm_1990)

Lindholm_2000 <- Lindholm %>% filter(Year == 2000) %>% 
  dplyr::select(c(1,3:68))
Lindholm_2000 <- right.matrix(Lindholm_2000)

Lindholm_2010 <- Lindholm %>% filter(Year == 2010) %>% 
  dplyr::select(c(1,3:68))
Lindholm_2010 <- right.matrix(Lindholm_2010)

Lindholm_met <- list(Lindholm_1940 = Lindholm_1940, Lindholm_1970 = Lindholm_1970,
                     Lindholm_1990 = Lindholm_1990, Lindholm_2000 = Lindholm_2000, 
                     Lindholm_2010 = Lindholm_2010)


Lindholm_env_1940 <- Lindholm_env %>% filter(time == "1940s") %>% 
  dplyr::select(c(1,5:12))

Lindholm_env_1940 <- right.matrix(Lindholm_env_1940)

Lindholm_env_1970 <- Lindholm_env %>% filter(time == "1970s") %>% 
  dplyr::select(c(1,5:12))

Lindholm_env_1970 <- right.matrix(Lindholm_env_1970)

Lindholm_env_1990 <- Lindholm_env %>% filter(time == "1990s") %>% 
  dplyr::select(c(1,5:12))

Lindholm_env_1990 <- right.matrix(Lindholm_env_1990)

Lindholm_env_2000 <- Lindholm_env %>% filter(time == "2000s") %>% 
  dplyr::select(c(1,5:12))

Lindholm_env_2000 <- right.matrix(Lindholm_env_2000)

Lindholm_env_2010 <- Lindholm_env %>% filter(time == "2010s") %>% 
  dplyr::select(c(1,5:12))

Lindholm_env_2010 <- right.matrix(Lindholm_env_2010)



Lindholm_env_ls <- list(Lindholm_env_1940 = Lindholm_env_1940, Lindholm_env_1970 = Lindholm_env_1970,
                        Lindholm_env_1990 = Lindholm_env_1990, Lindholm_env_2000 = Lindholm_env_2000, 
                        Lindholm_env_2010 = Lindholm_env_2010)





# network level analysis
sapply(Lindholm_met,function(x)networklevel(empty(x), index = c("connectance", "NODF")))

#               Lindholm_1940 Lindholm_1970 Lindholm_1990 Lindholm_2000 Lindholm_2010
# connectance     0.3557984     0.4020311     0.4074074      0.396412     0.3345091
# NODF           70.6633405    73.5329440    73.9240563     71.060736    64.7760313

sapply(Lindholm_met,function(x)computeModules(empty(x))@likelihood)

# Lindholm_1940 Lindholm_1970 Lindholm_1990 Lindholm_2000 Lindholm_2010 
# 0.2127745     0.1912877     0.1835688     0.2032776     0.2465945

sapply(Lindholm_met,function(x)sum(empty(x))) # links
# Lindholm_1940 Lindholm_1970 Lindholm_1990 Lindholm_2000 Lindholm_2010 
# 586           673           693           685           569



sapply(Lindholm_met,function(x)grouplevel(empty(x), index = c("number of species")))

sapply(Lindholm_met,function(x)specieslevel(empty(x), index = "degree", level = "lower"))
sapply(sapply(Lindholm_met,function(x)specieslevel(empty(x), index = "degree", level = "lower")),
       function(x)mean(x))


netstr_Lindholm <- as.data.frame(sapply(Lindholm_met,function(x)networklevel(empty(x), index = c("connectance", "NODF"))) %>% 
                                   rbind(.,sapply(Lindholm_met,function(x)computeModules(empty(x))@likelihood)) %>% 
                                   rbind(.,sapply(Lindholm_met,function(x)sum(empty(x)))) %>% 
                                   rbind(.,sapply(Lindholm_met,function(x)grouplevel(empty(x), 
                                                                                     index = "number of species", level = "higher"))) %>% 
                                   rbind(.,sapply(Lindholm_met,function(x)mean(rowSums(empty(x))))))

str(netstr_Lindholm)

netstr_Lindholm$Index <-c("connectance", "NODF", "Modularity", "number of links", 
                          "gamma species richness", "alpha species richness")

netstr_Lindholm_plot <- reshape2::melt(netstr_Lindholm,id = "Index") %>% 
  ggplot(.,aes(x = variable, y = value, color = Index, group = Index)) +
  geom_line(size = 1)+geom_point(size = 3) + 
  scale_y_log10() + scale_x_discrete(breaks = c("Lindholm_1940", "Lindholm_1970",
                                                "Lindholm_1990", "Lindholm_2000", 
                                                "Lindholm_2010"), 
                                     labels = c("1940s", "1970s", "1990s",
                                                "2000s", "2010s")) +
  scale_color_brewer(palette = "Dark2") + 
  xlab("") + ylab("Species richness and network structure index") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = c(0.35, 0.35)) + 
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))

netstr_Lindholm_plot


# node level analysis

degree_Lindholm_time <- specieslevel(empty(Lindholm_1940), index = "degree", level = "higher") %>% 
  mutate(node = rownames(.), time = "1940s", type = "plant") %>% 
  rbind(.,specieslevel(empty(Lindholm_1970), index = "degree", level = "higher") %>% 
          mutate(node = rownames(.), time = "1970s", type = "plant")) %>% 
  rbind(.,specieslevel(empty(Lindholm_1990), index = "degree", level = "higher") %>% 
          mutate(node = rownames(.), time = "1990s", type = "plant")) %>% 
  rbind(.,specieslevel(empty(Lindholm_2000), index = "degree", level = "higher") %>% 
          mutate(node = rownames(.), time = "2000s", type = "plant")) %>% 
  rbind(.,specieslevel(empty(Lindholm_2010), index = "degree", level = "higher") %>% 
          mutate(node = rownames(.), time = "2010s", type = "plant")) %>% 
  rbind(.,specieslevel(empty(Lindholm_1940), index = "degree", level = "lower") %>% 
          mutate(node = rownames(.), time = "1940s", type = "lake")) %>%
  rbind(.,specieslevel(empty(Lindholm_1970), index = "degree", level = "lower") %>% 
          mutate(node = rownames(.), time = "1970s", type = "lake")) %>%
  rbind(.,specieslevel(empty(Lindholm_1990), index = "degree", level = "lower") %>% 
          mutate(node = rownames(.), time = "1990s", type = "lake")) %>%
  rbind(.,specieslevel(empty(Lindholm_2000), index = "degree", level = "lower") %>% 
          mutate(node = rownames(.), time = "2000s", type = "lake")) %>%
  rbind(.,specieslevel(empty(Lindholm_2010), index = "degree", level = "lower") %>% 
          mutate(node = rownames(.), time = "2010s", type = "lake"))


degree_Lindholm_time$time <- as.factor(degree_Lindholm_time$time)

fit_degree <- aov(log(degree) ~ time, data = degree_Lindholm_time %>% 
                    filter(type == "lake"))

library(multcomp)

tuk_degree <- glht(fit_degree, linfct = mcp(time = "Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk_degree, level=.05),col="lightgrey")
par(op)


fit2_degree <- aov(log(degree) ~ time, data = degree_Lindholm_time %>% 
                     filter(type == "plant"))


tuk2_degree <- glht(fit2_degree, linfct = mcp(time = "Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk2_degree, level = .05),col = "lightgrey")
par(op)

0.825

degree_time_boxplot <- ggplot(degree_Lindholm_time,aes(x = time, y = degree,fill = type))+
  geom_boxplot()+
  annotate("text", x = 0.825, y = 45, label = "a",
           col = "black", size = 6) +
  annotate("text", x = 1.825, y = 45, label = "a",
           col = "black", size = 6) +
  annotate("text", x = 2.825, y = 45, label = "a",
           col = "black", size = 6) +
  annotate("text", x = 3.825, y = 45, label = "a",
           col = "black", size = 6) +
  annotate("text", x = 4.825, y = 45, label = "a",
           col = "black", size = 6) +
  annotate("text", x = 0.825+0.375, y = 28, label = "A",
           col = "black", size = 6)+
  annotate("text", x = 1.825+0.375, y = 28, label = "A",
           col = "black", size = 6)+
  annotate("text", x = 2.825+0.375, y = 28, label = "A",
           col = "black", size = 6)+
  annotate("text", x = 3.825+0.375, y = 28, label = "A",
           col = "black", size = 6)+
  annotate("text", x = 4.825+0.375, y = 28, label = "A",
           col = "black", size = 6)+
  ylab("Degree") + 
  xlab("") + 
  scale_fill_manual(values = c("#00A087B2", "#8491B4B2"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))

degree_time_boxplot




degree_sp_Lindholm <- specieslevel(empty(Lindholm_1940), index = "degree", level = "higher") %>% 
  mutate(species = rownames(.), time = "1940s") %>% 
  rbind(.,specieslevel(empty(Lindholm_1970), index = "degree", level = "higher") %>% 
          mutate(species = rownames(.), time = "1970s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_1990), index = "degree", level = "higher") %>% 
          mutate(species = rownames(.), time = "1990s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_2000), index = "degree", level = "higher") %>% 
          mutate(species = rownames(.), time = "2000s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_2010), index = "degree", level = "higher") %>% 
          mutate(species = rownames(.), time = "2010s")) %>% 
  group_by(species) %>% 
  summarise(mu = mean(degree),
            sd = sd(degree),
            n = n()) %>% 
  mutate(cv = sd/mu)

cv_degree_density_plot <- degree_Lindholm_time %>% 
  group_by(node,type) %>% 
  summarise(mu = mean(degree),
            sd = sd(degree),
            n = n()) %>% 
  mutate(cv = sd/mu) %>% 
  ggplot(., aes(x = cv))+
  geom_density(aes(fill=type), color = NA,alpha = 0.5)+
  ylab("Density") + 
  xlab("Coefficient of variation of degree of each node") + 
  scale_fill_manual(values = c("#00A087B2", "#8491B4B2"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))

cv_degree_density_plot


degree_patch_Lindholm_time <- specieslevel(empty(Lindholm_1940), index = "degree", level = "lower") %>% 
  mutate(patch = rownames(.), time = "1940s") %>% 
  rbind(.,specieslevel(empty(Lindholm_1970), index = "degree", level = "lower") %>% 
          mutate(patch = rownames(.), time = "1970s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_1990), index = "degree", level = "lower") %>% 
          mutate(patch = rownames(.), time = "1990s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_2000), index = "degree", level = "lower") %>% 
          mutate(patch = rownames(.), time = "2000s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_2010), index = "degree", level = "lower") %>% 
          mutate(patch = rownames(.), time = "2010s")) %>%
  full_join(.,Lindholm_env, by = c("patch", "time"))



degree_patch_Lindholm <- specieslevel(empty(Lindholm_1940), index = "degree", level = "lower") %>% 
  mutate(patch = rownames(.), time = "1940s") %>% 
  rbind(.,specieslevel(empty(Lindholm_1970), index = "degree", level = "lower") %>% 
          mutate(patch = rownames(.), time = "1970s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_1990), index = "degree", level = "lower") %>% 
          mutate(patch = rownames(.), time = "1990s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_2000), index = "degree", level = "lower") %>% 
          mutate(patch = rownames(.), time = "2000s")) %>% 
  rbind(.,specieslevel(empty(Lindholm_2010), index = "degree", level = "lower") %>% 
          mutate(patch = rownames(.), time = "2010s")) %>% 
  group_by(patch) %>% 
  summarise(mu = mean(degree),
            sd = sd(degree),
            n = n()) %>% 
  mutate(cv = sd/mu)


cnodf_1940 <- nestedcontribution(empty(Lindholm_1940), nsimul = 999)
cnodf_1970 <- nestedcontribution(empty(Lindholm_1970), nsimul = 999)
cnodf_1990 <- nestedcontribution(empty(Lindholm_1990), nsimul = 999)
cnodf_2000 <- nestedcontribution(empty(Lindholm_2000), nsimul = 999)
cnodf_2010 <- nestedcontribution(empty(Lindholm_2010), nsimul = 999)

save(cnodf_1940, cnodf_1970, cnodf_1990,cnodf_2000,cnodf_2010,
     file = "C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Data/cnodf_lindholm.Rdata")


cnodf_sp_Lindholm <- cnodf_1940$`higher level` %>% 
  mutate(species = rownames(.), time = "1940s") %>% 
  rbind(.,cnodf_1970$`higher level` %>% 
          mutate(species = rownames(.), time = "1970s")) %>% 
  rbind(.,cnodf_1990$`higher level` %>% 
          mutate(species = rownames(.), time = "1990s")) %>% 
  rbind(.,cnodf_2000$`higher level` %>% 
          mutate(species = rownames(.), time = "2000s")) %>% 
  rbind(.,cnodf_2010$`higher level` %>% 
          mutate(species = rownames(.), time = "2010s")) %>% 
  group_by(species) %>% 
  summarise(mu = mean(nestedcontribution),
            sd = sd(nestedcontribution),
            n = n()) %>% 
  mutate(cv = sd/mu)

hist(cnodf_sp_Lindholm$cv, freq = F)

plot(density(na.omit(cnodf_sp_Lindholm$cv)))

hist(degree_sp_Lindholm$cv, freq = F)

ggplot(cnodf_sp_Lindholm, aes(x = reorder(species, -cv), y = cv))+
  geom_bar(stat = "identity")

plot(degree_sp_Lindholm$cv)

plot(cnodf_sp_Lindholm$cv)
abline(h = -1)
abline(h = 0)
abline(h = 1)


cnodf_patch_Lindholm <- cnodf_1940$`lower level` %>% 
  mutate(species = rownames(.), time = "1940s") %>% 
  rbind(.,cnodf_1970$`lower level` %>% 
          mutate(species = rownames(.), time = "1970s")) %>% 
  rbind(.,cnodf_1990$`lower level` %>% 
          mutate(species = rownames(.), time = "1990s")) %>% 
  rbind(.,cnodf_2000$`lower level` %>% 
          mutate(species = rownames(.), time = "2000s")) %>% 
  rbind(.,cnodf_2010$`lower level` %>% 
          mutate(species = rownames(.), time = "2010s")) %>% 
  group_by(species) %>% 
  summarise(mu = mean(nestedcontribution),
            sd = sd(nestedcontribution),
            n = n()) %>% 
  mutate(cv = sd/mu)


cnodf_patch_Lindholm_time <- cnodf_1940$`lower level` %>% 
  mutate(patch = rownames(.), time = "1940s") %>% 
  rbind(., cnodf_1970$`lower level` %>% 
          mutate(patch = rownames(.), time = "1970s")) %>% 
  rbind(., cnodf_1990$`lower level` %>% 
          mutate(patch = rownames(.), time = "1990s")) %>% 
  rbind(., cnodf_2000$`lower level` %>% 
          mutate(patch = rownames(.), time = "2000s")) %>% 
  rbind(., cnodf_2010$`lower level` %>% 
          mutate(patch = rownames(.), time = "2010s")) %>%
  full_join(., Lindholm_env, by = c("patch", "time"))




library(ggpubr)
ggplot(cnodf_patch_Lindholm_time, aes(x = pH, y = nestedcontribution))+
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y = 4, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 3.8, aes(label = ..adj.rr.label..)) +
  facet_grid(~ time)

names(cnodf_patch_Lindholm_time)

rf_cnodf2



cnodf_Lindholm_time <- cnodf_1940$`lower level` %>% 
  mutate(node = rownames(.), time = "1940s", type = "lake") %>% 
  rbind(.,cnodf_1970$`lower level` %>% 
          mutate(node = rownames(.), time = "1970s", type = "lake")) %>% 
  rbind(.,cnodf_1990$`lower level` %>% 
          mutate(node = rownames(.), time = "1990s", type = "lake")) %>% 
  rbind(.,cnodf_2000$`lower level` %>% 
          mutate(node = rownames(.), time = "2000s", type = "lake")) %>% 
  rbind(.,cnodf_2010$`lower level` %>% 
          mutate(node = rownames(.), time = "2010s", type = "lake")) %>%
  rbind(.,cnodf_1940$`higher level` %>% 
          mutate(node = rownames(.), time = "1940s", type = "plant")) %>%
  rbind(.,cnodf_1970$`higher level` %>% 
          mutate(node = rownames(.), time = "1970s", type = "plant")) %>% 
  rbind(.,cnodf_1990$`higher level` %>% 
          mutate(node = rownames(.), time = "1990s", type = "plant")) %>% 
  rbind(.,cnodf_2000$`higher level` %>% 
          mutate(node = rownames(.), time = "2000s", type = "plant")) %>% 
  rbind(.,cnodf_2010$`higher level` %>% 
          mutate(node = rownames(.), time = "2010s", type = "plant"))



summary(aov(nestedcontribution ~ time, data = cnodf_Lindholm_time %>% 
              filter(type == "plant")))

cnodf_Lindholm_time$time <- as.factor(cnodf_Lindholm_time$time)

fit <- aov(nestedcontribution ~ time, data = cnodf_Lindholm_time %>% 
             filter(type == "lake"))
library(multcomp)

tuk <- glht(fit, linfct = mcp(time = "Tukey"))
op <- par(mar=c(5, 4 ,6, 2))
plot(cld(tuk, level = .05), col = "lightgrey")
par(op)


fit2 <- aov(nestedcontribution ~ time, data = cnodf_Lindholm_time %>% 
              filter(type == "plant"))


tuk2 <- glht(fit2, linfct = mcp(time = "Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk2, level=.05), col = "lightgrey")
par(op)



cnodf_time_boxplot <- ggplot(cnodf_Lindholm_time, aes(x = time, y = nestedcontribution,fill = type))+
  geom_boxplot()+
  annotate("text", x = 0.825, y = 3.5, label = "ab",
           col = "black", size = 6) +
  annotate("text", x = 1.825, y = 3.5, label = "ab",
           col = "black", size = 6) +
  annotate("text", x = 2.825, y = 3.5, label = "b",
           col = "black", size = 6) +
  annotate("text", x = 3.825, y = 3.5, label = "a",
           col = "black", size = 6) +
  annotate("text", x = 4.825, y = 3.5, label = "a",
           col = "black", size = 6) +
  annotate("text", x = 0.825 + 0.375, y = 3.2, label = "A",
           col = "black", size = 6) +
  annotate("text", x = 1.825 + 0.375, y = 3.2, label = "A",
           col = "black", size = 6) +
  annotate("text", x = 2.825 + 0.375, y = 3.2, label = "A",
           col = "black", size = 6) +
  annotate("text", x = 3.825 + 0.375, y = 3.2, label = "A",
           col = "black", size = 6) +
  annotate("text", x = 4.825 + 0.375, y = 3.2, label = "A",
           col = "black", size = 6) +
  ylab("Contribution to nestedness (cnodf)") + 
  xlab("") + 
  scale_fill_manual(values = c("#00A087B2", "#8491B4B2"))+
  theme_bw() + 
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))

cnodf_time_boxplot



## Adding the minimum to everyone so the cv makes sense, cf. Gaiarsa et al. 2021 NEE

smallest_cnodf <- abs(min(cnodf_Lindholm_time$nestedcontribution))
cnodf_Lindholm_time$nestedcontribution2 <- cnodf_Lindholm_time$nestedcontribution + smallest_cnodf


cv_cnodf_Lind <- cnodf_Lindholm_time %>% 
  group_by(node,type) %>% 
  summarise(mu = mean(nestedcontribution2),
            sd = sd(nestedcontribution2),
            n = n()) %>% 
  mutate(cv = sd/mu,Index = "cnodf") %>% 
  filter(n > 2)

cv_cnodf_Lind$type <- as.factor(cv_cnodf_Lind$type)

cv_degree_Lind <- degree_Lindholm_time %>% 
  group_by(node,type) %>% 
  summarise(mu = mean(degree),
            sd = sd(degree),
            n = n()) %>% 
  mutate(cv = sd/mu,Index = "degree") %>% 
  filter(n > 2)

summary(aov(cv ~ type, cv_cnodf_Lind))
summary(aov(cv ~ type, cv_degree_Lind))
cv_degree_Lind$type <- as.factor(cv_degree_Lind$type)


cv_all_Lind <- cnodf_Lindholm_time %>% 
  group_by(node, type) %>% 
  summarise(mu = mean(nestedcontribution2),
            sd = sd(nestedcontribution2),
            n = n()) %>% 
  mutate(cv = sd/mu, Index = "cnodf") %>% 
  rbind(., degree_Lindholm_time %>% 
          group_by(node, type) %>% 
          summarise(mu = mean(degree),
                    sd = sd(degree),
                    n = n()) %>% 
          mutate(cv = sd/mu, Index = "degree")) %>% 
  filter(n > 2)

cv_all_Lind$variable <- paste(cv_all_Lind$type, cv_all_Lind$Index, sep = "_")

cv_all_Lind$variable <- as.factor(cv_all_Lind$variable)

fit3 <- aov(cv ~ variable, data = cv_all_Lind)

summary(fit3)


fit3 <- aov(cv~type, cv_cnodf_Lind)

tuk3 <- glht(fit3, linfct = mcp(type = "Tukey"))
op <- par(mar=c(5,4,6,2))
plot(cld(tuk3, level=.05),col="lightgrey")
par(op)




cv_boxplot <- cnodf_Lindholm_time %>% 
  group_by(node,type) %>% 
  summarise(mu = mean(nestedcontribution2),
            sd = sd(nestedcontribution2),
            n = n()) %>% 
  mutate(cv = sd/mu,Index = "cnodf") %>% 
  rbind(.,degree_Lindholm_time %>% 
          group_by(node,type) %>% 
          summarise(mu = mean(degree),
                    sd = sd(degree),
                    n = n()) %>% 
          mutate(cv = sd/mu,Index = "degree")) %>% 
  filter(n > 2) %>% 
  ggplot(.,aes(x = Index, y = cv, fill = type)) +
  geom_boxplot()+
  xlab("") + 
  ylab("Structural variability") + 
  scale_fill_manual(values = c("#00A087B2", "#8491B4B2"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))

cv_boxplot


cv_boxplot2 <- cnodf_Lindholm_time %>% 
  group_by(node,type) %>% 
  summarise(mu = mean(nestedcontribution2),
            sd = sd(nestedcontribution2),
            n = n()) %>% 
  mutate(cv = sd/mu,Index = "cnodf") %>% 
  filter(n > 2) %>% 
  ggplot(.,aes(x = type, y = cv, fill = type)) +
  geom_boxplot(width = 0.2)+
  xlab("") + 
  ylab("Structural variability") + 
  scale_fill_manual(values = c("#00A087B2", "#8491B4B2"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))
cv_boxplot2


my_comparisons <- list(c("lake_cnodf", "plant_cnodf"),
                       c("lake_cnodf", "lake_degree"),
                       c("lake_cnodf", "plant_degree"),
                       c("lake_degree", "plant_degree"),
                       c("plant_degree", "plant_cnodf"))

ggplot(cv_all_Lind,aes(x = variable, y = cv, fill = type)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 0.8)                      



cv_cnodf_density_plot <- cnodf_Lindholm_time %>% 
  group_by(node,type) %>% 
  summarise(mu = mean(nestedcontribution),
            sd = sd(nestedcontribution),
            n = n()) %>% 
  mutate(cv = sd/mu) %>% 
  ggplot(., aes(x = cv)) +
  geom_density(aes(fill = type), color = NA,alpha = 0.5) +
  ylab("Density") + 
  xlab("CV of nestedness contribution of each node") + 
  scale_fill_manual(values = c("#00A087B2", "#8491B4B2")) +
  theme_bw() +
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))

cv_cnodf_density_plot




library(corrplot)
M = cor(cnodf_patch_Lindholm_time[,c(6:13)])
testRes = cor.mtest(cnodf_patch_Lindholm_time[,c(6:13)], conf.level = 0.95)

corrplot(M, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20')


library(MuMIn)

fit_lmm <- lmer(nestedcontribution ~ Elevation + Depth + Area + pH +
                  Secchi + Agriculture+ Built.area + Ditches+
                  (1|time), data = cnodf_patch_Lindholm_time,
                REML = FALSE,
                na.action = "na.fail")

summary(fit_lmm)
vif(fit_lmm)

# Elevation       Depth        Area          pH      Secchi 
# 8.204269    2.254344    1.694982    2.596196    3.004046 
# Agriculture  Built.area     Ditches 
# 4.558729    2.609783    1.300597

# removing elevation
fit_lmm_update <- lmer(nestedcontribution ~ Depth + Area + pH +
                         Secchi + Agriculture+ Built.area + Ditches+
                         (1|time), data = cnodf_patch_Lindholm_time,
                       REML = FALSE,
                       na.action = "na.fail")

summary(fit_lmm_update)
vif(fit_lmm_update)


plot(fit_lmm_update)

car::Anova(fit_lmm_update)

plot(DHARMa::simulateResiduals(fit_lmm_update, plot = FALSE)) # model check, 

res <- dredge(fit_lmm_update)

top_mod <- subset(res, delta <= 2, recalc.weights = FALSE)
top_mod

summary(model.avg(res, revised.var = FALSE))

MuMIn::importance(res)
importance(top_mod)

visreg(model.avg(top_mod, revised.var = FALSE))



fit_lmm_1 <- lmer(nestedcontribution ~ pH +
                    (1|time), data = cnodf_patch_Lindholm_time,
                  na.action = "na.fail")


car::Anova(fit_lmm_1)
plot(fit_lmm_1)
r.squaredGLMM(fit_lmm_1)

plot(DHARMa::simulateResiduals(fit_lmm_1, plot = FALSE))

fit_lmm_2 <- lmer(nestedcontribution ~ pH + Secchi + 
                    (1|time), data = cnodf_patch_Lindholm_time,
                  na.action = "na.fail")

car::Anova(fit_lmm_2)
plot(DHARMa::simulateResiduals(fit_lmm_2, plot = FALSE))

fit_lmm_3 <- lmer(nestedcontribution ~ Built.area +  pH + 
                    (1|time), data = cnodf_patch_Lindholm_time,
                  na.action = "na.fail")
summary(fit_lmm_3)

car::Anova(fit_lmm_3)


fit_lmm_4 <- lmer(nestedcontribution ~ Elevation + 
                    (1|time), data = cnodf_patch_Lindholm_time,
                  na.action = "na.fail")
summary(fit_lmm_4)

car::Anova(fit_lmm_4)


ggplot() +
  geom_point(data = cnodf_patch_Lindholm_time,
             aes(x = pH, y = nestedcontribution, color = time),size = 4)+
  geom_smooth(data = cnodf_patch_Lindholm_time,
              aes(x = pH, y = nestedcontribution, color = time),method = "lm", se = F) +
  geom_ribbon(data = predict_fun(model = fit_lmm_1, 
                                 variable = "pH", 
                                 fit = "nestedcontribution"), 
              aes(x = pH, ymin = lower, ymax = upper),
              color = NA, alpha=.5) +
  geom_line(data = predict_fun(model = fit_lmm_1, 
                               variable = "pH", 
                               fit = "nestedcontribution"), 
            aes(x = pH, y = nestedcontribution), colour="#1C86EE",
            size=1, linetype = 1)+
  scale_color_brewer(palette = "Dark2") +
  xlab("pH") + ylab("Contribution to nestedness") + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))

ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/lindholm_cnodf_pH_figure.png",
       height = 6, width = 6.6, units = "in", dpi = 600)
ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/lindholm_cnodf_pH_figure.pdf",
       height = 6, width = 6.6, units = "in")


ggplot() +
  geom_point(data = cnodf_patch_Lindholm_time,
             aes(x = Elevation, y = nestedcontribution, color = time),size = 4)+
  geom_smooth(data = cnodf_patch_Lindholm_time,
              aes(x = Elevation, y = nestedcontribution, color = time),method = "lm", se = F) +
  geom_ribbon(data = predict_fun(model = fit_lmm_4, 
                                 variable = "Elevation", 
                                 fit = "nestedcontribution"), 
              aes(x = Elevation, ymin = lower, ymax = upper),
              color = NA, alpha=.5)+
  geom_line(data = predict_fun(model = fit_lmm_4, 
                               variable = "Elevation", 
                               fit = "nestedcontribution"), 
            aes(x = Elevation, y = nestedcontribution), colour="#1C86EE",
            size = 1, linetype = 1)+
  scale_color_brewer(palette = "Dark2") +
  xlab("Elevation") + ylab("Contribution to nestedness") + 
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"))


ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/lindholm_cnodf_elev_figure.png",
       height = 6, width = 6.6, units = "in", dpi = 600)

ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/lindholm_cnodf_elev_figure.pdf",
       height = 6, width = 6.6, units = "in")




fit_glmm <- glmer(degree ~ Elevation+Depth + Area+pH+
                    Secchi + Agriculture+ Built.area+ Ditches+
                    (1 | time), data = degree_patch_Lindholm_time,
                  family = poisson(),
                  na.action = "na.fail")
summary(fit_glmm)

visreg::visreg(fit_glmm, scale = "response", type = "contrast")

plot(fit_glmm)
car::vif(fit_glmm)

car::Anova(fit_lmm)

plot(DHARMa::simulateResiduals(fit_glmm, plot = FALSE)) # model check, 


specieslevel(empty(Lindholm_1990), index = "degree", level = "higher")
specieslevel(empty(Lindholm_2000), index = "degree", level = "higher")
specieslevel(empty(Lindholm_2010), index = "degree", level = "higher")


Lindholm_beta_f_l <- list()
Lindholm_beta_0_l <- list()
#env_f_1 <- list()
#env_0_1 <- list()


for (i in 1:4) {
  beta <- beta_temporal_metacomm (Lindholm_met[[i]], Lindholm_met[[i+1]])
  beta$time <- paste(names(Lindholm_met)[i],"<->", names(Lindholm_met)[i + 1],sep = "")
  beta$Elevation <- dist(rbind(Lindholm_env_ls[[i]][,"Elevation"], Lindholm_env_ls[[i + 1]][, "Elevation"]))
  beta$Depth <- dist(rbind(Lindholm_env_ls[[i]][, "Depth"], Lindholm_env_ls[[i + 1]][, "Depth"]))
  beta$Area <- dist(rbind(Lindholm_env_ls[[i]][, "Area"], Lindholm_env_ls[[i + 1]][,"Area"]))
  beta$pH <- dist(rbind(Lindholm_env_ls[[i]][, "pH"],Lindholm_env_ls[[i + 1]][,"pH"]))
  beta$Secchi <- dist(rbind(Lindholm_env_ls[[i]][, "Secchi"],Lindholm_env_ls[[i + 1]][,"Secchi"]))
  beta$Agriculture <- dist(rbind(Lindholm_env_ls[[i]][, "Agriculture"],Lindholm_env_ls[[i + 1]][,"Agriculture"]))
  beta$Built.area <- dist(rbind(Lindholm_env_ls[[i]][, "Built.area"],Lindholm_env_ls[[i + 1]][,"Built.area"]))
  beta$Ditches <- dist(rbind(Lindholm_env_ls[[i]][, "Ditches"],Lindholm_env_ls[[i + 1]][,"Ditches"]))
  
  beta0 <- beta_temporal_metacomm (Lindholm_met[[1]], Lindholm_met[[i + 1]])
  beta0$time <- paste(names(Lindholm_met)[1], "<->", names(Lindholm_met)[i + 1],sep = "")
  beta0$Elevation <- dist(rbind(Lindholm_env_ls[[1]][, "Elevation"], Lindholm_env_ls[[i + 1]][, "Elevation"]))
  beta0$Depth <- dist(rbind(Lindholm_env_ls[[1]][,"Depth"], Lindholm_env_ls[[i + 1]][, "Depth"]))
  beta0$Area <- dist(rbind(Lindholm_env_ls[[1]][,"Area"], Lindholm_env_ls[[i + 1]][, "Area"]))
  beta0$pH <- dist(rbind(Lindholm_env_ls[[1]][,"pH"], Lindholm_env_ls[[i + 1]][, "pH"]))
  beta0$Secchi <- dist(rbind(Lindholm_env_ls[[1]][, "Secchi"], Lindholm_env_ls[[i + 1]][, "Secchi"]))
  beta0$Agriculture <- dist(rbind(Lindholm_env_ls[[1]][, "Agriculture"], Lindholm_env_ls[[i + 1]][, "Agriculture"]))
  beta0$Built.area <- dist(rbind(Lindholm_env_ls[[1]][, "Built.area"], Lindholm_env_ls[[i + 1]][, "Built.area"]))
  beta0$Ditches <- dist(rbind(Lindholm_env_ls[[1]][, "Ditches"], Lindholm_env_ls[[i + 1]][, "Ditches"]))
  
  
  Lindholm_beta_f_l[[i]] <- beta
  Lindholm_beta_0_l[[i]] <- beta0
  
}




Lindholm_beta_f <- do.call("rbind", Lindholm_beta_f_l)
Lindholm_beta_f
write.csv(Lindholm_beta_f,"C:/Dataset/Temporal_metacommunity/Res/Lindholm_beta_f.csv")


Lindholm_beta_0 <- do.call("rbind", Lindholm_beta_0_l)
Lindholm_beta_0
write.csv(Lindholm_beta_0,"C:/Dataset/Temporal_metacommunity/Res/Lindholm_beta_0.csv")


Lindholm_beta_f_l <- list()
Lindholm_beta_0_l <- list()
#env_f_1 <- list()
#env_0_1 <- list()


for (i in 1:4) {
  beta <- beta_temporal_metacomm (Lindholm_met[[i]], Lindholm_met[[i + 1]])
  beta$time <- paste(names(Lindholm_met)[i], "<->", names(Lindholm_met)[i + 1], sep = "")
  beta$Elevation <- mean(Lindholm_env_ls[[i + 1]][, "Elevation"] - Lindholm_env_ls[[i]][, "Elevation"])
  beta$Depth <- mean(Lindholm_env_ls[[i + 1]][, "Depth"] - Lindholm_env_ls[[i]][, "Depth"])
  beta$Area <-mean(Lindholm_env_ls[[i + 1]][, "Area"] - Lindholm_env_ls[[i]][, "Area"])
  beta$pH <- mean(Lindholm_env_ls[[i + 1]][, "pH"] - Lindholm_env_ls[[i]][, "pH"])
  beta$Secchi <- mean(Lindholm_env_ls[[i + 1]][, "Secchi"] - Lindholm_env_ls[[i]][, "Secchi"])
  beta$Agriculture <- mean(Lindholm_env_ls[[i + 1]][, "Agriculture"] - Lindholm_env_ls[[i]][, "Agriculture"])
  beta$Built.area <-mean(Lindholm_env_ls[[i + 1]][, "Built.area"] - Lindholm_env_ls[[i]][, "Built.area"])
  beta$Ditches <- mean(Lindholm_env_ls[[i + 1]][, "Ditches"] - Lindholm_env_ls[[i]][, "Ditches"])
  
  beta0 <- beta_temporal_metacomm (Lindholm_met[[1]], Lindholm_met[[i + 1]])
  beta0$time <- paste(names(Lindholm_met)[1], "<->", names(Lindholm_met)[i + 1], sep = "")
  
  beta0$Elevation <- mean(Lindholm_env_ls[[i + 1]][, "Elevation"] - Lindholm_env_ls[[1]][, "Elevation"])
  beta0$Depth <- mean(Lindholm_env_ls[[i + 1]][, "Depth"] - Lindholm_env_ls[[1]][, "Depth"])
  beta0$Area <- mean(Lindholm_env_ls[[i + 1]][, "Area"] - Lindholm_env_ls[[1]][, "Area"])
  beta0$pH <- mean(Lindholm_env_ls[[i + 1]][, "pH"] - Lindholm_env_ls[[1]][, "pH"])
  beta0$Secchi <- mean(Lindholm_env_ls[[i + 1]][, "Secchi"] - Lindholm_env_ls[[1]][, "Secchi"])
  beta0$Agriculture <- mean(Lindholm_env_ls[[i + 1]][, "Agriculture"] - Lindholm_env_ls[[1]][, "Agriculture"])
  beta0$Built.area <- mean(Lindholm_env_ls[[i + 1]][, "Built.area"] - Lindholm_env_ls[[1]][, "Built.area"])
  beta0$Ditches <- mean(Lindholm_env_ls[[i + 1]][, "Ditches"] - Lindholm_env_ls[[1]][, "Ditches"])
  
  
  Lindholm_beta_f_l[[i]] <- beta
  Lindholm_beta_0_l[[i]] <- beta0
  
}




Lindholm_beta_f <- do.call("rbind", Lindholm_beta_f_l)
Lindholm_beta_f
write.csv(Lindholm_beta_f,"C:/Dataset/Temporal_metacommunity/Res/Lindholm_beta_f.csv")


Lindholm_beta_0 <- do.call("rbind", Lindholm_beta_0_l)
Lindholm_beta_0
write.csv(Lindholm_beta_0,"C:/Dataset/Temporal_metacommunity/Res/Lindholm_beta_0.csv")




ggplot(Lindholm_beta_f %>% 
         filter(Component != "beta_Landscape" & Component != "beta_RL") %>% 
         reshape2::melt(.,id = c("Component", "beta", "time", "Elevation")),
       aes(x = value, y = beta)) +
  geom_point(aes(color = factor(time,levels =  c("Lindholm_1940<->Lindholm_1970",
                                                 "Lindholm_1970<->Lindholm_1990",
                                                 "Lindholm_1990<->Lindholm_2000",
                                                 "Lindholm_2000<->Lindholm_2010"),
                                labels = c("1940-1970",
                                           "1970-1990",
                                           "1990-2000",
                                           "2000-2010"))), size = 4) +
  geom_smooth(method = "lm", se = F, color = "grey50", size = 1.5)+
  geom_path(arrow = arrow(ends = "last",length = unit(0.1,"inches")), size = 1, color = "blue")+
  xlab("△ value") + ylab("β-diversity") + 
  facet_grid(Component ~ variable, scales = "free")+
  theme_bw() + 
  theme(legend.title = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, family = "serif"),
        strip.text = element_text(size = 12, family = "serif"))

ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/lindholm_beta_env_figure.pdf",
       height = 7*1.5, width = 10*1.5, units = "in")


Lindholm_beta_f$Depth

ggplot(Lindholm_beta_0 %>% 
         filter(Component != "beta_Landscape" & Component != "beta_RL") %>% 
         reshape2::melt(.,id = c("Component", "beta", "time", "Elevation")),
       aes(x = value, y = beta)) +
  geom_point(aes(color = factor(time, levels = c("Lindholm_1940<->Lindholm_1970",
                                                 "Lindholm_1940<->Lindholm_1990",
                                                 "Lindholm_1940<->Lindholm_2000",
                                                 "Lindholm_1940<->Lindholm_2010"),
                                labels = c("1940-1970",
                                           "1940-1990",
                                           "1940-2000",
                                           "1940-2010")))) +
  geom_smooth(method = "lm")+
  facet_grid(Component ~ variable, scales = "free")+
  theme_bw() + 
  theme(legend.title = element_blank())




lindholm_plot_a <- ggplot(Lindholm_beta_f,aes(x = time, y = beta, color = Component, group = Component))+
  geom_point(size = 2) + geom_line()+
  scale_color_brewer(palette = "Dark2")+
  scale_x_discrete(breaks = c("Lindholm_1940<->Lindholm_1970",
                              "Lindholm_1970<->Lindholm_1990",
                              "Lindholm_1990<->Lindholm_2000",
                              "Lindholm_2000<->Lindholm_2010"),
                   labels = c("1940-1970",
                              "1970-1990",
                              "1990-2000",
                              "2000-2010"))+
  ylab("β-diversity") + xlab("")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 12, family = "serif"))+
  theme(legend.position = "none")

lindholm_plot_a


lindholm_plot_b <- ggplot(Lindholm_beta_0,aes(x = time, y = beta, color = Component, group = Component))+
  geom_point(size = 2) + geom_line()+
  scale_color_brewer(palette = "Dark2")+
  scale_x_discrete(breaks = c("Lindholm_1940<->Lindholm_1970",
                              "Lindholm_1940<->Lindholm_1990",
                              "Lindholm_1940<->Lindholm_2000",
                              "Lindholm_1940<->Lindholm_2010"),
                   labels = c("1940-1970",
                              "1940-1990",
                              "1940-2000",
                              "1940-2010"))+
  ylab("β-diversity") + xlab("")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.title.y = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 12, family = "serif"))

lindholm_plot_b

lindholm_plot <- lindholm_plot_a + lindholm_plot_b + patchwork::plot_layout(ncol = 2)
lindholm_plot

ggsave("C:/Dataset/Temporal_metacommunity/Figs/lindholm_plot.pdf",
       plot = lindholm_plot,
       height = 4.31*2.5, width = 11.39*2.5, units = "cm")

ggsave("C:/Dataset/Temporal_metacommunity/Figs/lindholm_plot.png",
       plot = lindholm_plot,
       height = 4.31*2.5, width = 11.39*2.5, units = "cm")



netstr_Lindholm_plot + degree_time_boxplot +
  cnodf_time_boxplot + cv_boxplot2 + lindholm_plot_a + lindholm_plot_b +
  plot_layout(ncol = 2)

ggsave("C:/Dataset/Temporal_metacommunity/Temporal_beta_metacommunity/Figs/lindholm_figure_v2.pdf",
       height = 7*1.5, width = 8*1.5, units = "in")


## network----



ig <- list()
for (i in 1:5) {
  mt <- Lindholm_met[[i]]
  igt <- graph_from_incidence_matrix(mt, weighted = TRUE)
  
  vertex_attr(igt)$color <- rep("#8491B4B2", length(V(igt)))
  vertex_attr(igt)$color[grep(pattern = "FALSE", vertex_attr(igt)$type)] <- "#00A087B2"
    
  vertex_attr(igt)$color
  vertex_attr(igt)$shape<-rep("circle", length(V(igt)))
  vertex_attr(igt)$shape[grep(pattern = "FALSE", vertex_attr(igt)$type)]<- "square"
  ig[[i]] <- igt
}

l <- layout_nicely(ig[[1]])

years <- c(1940, 1970, 1990, 2000, 2010)


for(i in 1:5) {
  plot(
    ig[[i]],
    vertex.color = vertex_attr(ig[[i]])$cor,
    vertex.label = NA,
    vertex.size = 8,
    edge.width = 1,
    edge.color = "grey50",
    edge.curved = 0.3,
    layout = l
  )
  
  pdf(paste('C:/Dataset/Temporal_metacommunity/Figs/network_',years[i],'.pdf',sep=''),width = 6,height = 6)
  plot(
    ig[[i]],
    vertex.color = vertex_attr(ig[[i]])$cor,
    vertex.label = NA,
    vertex.size = 8,
    edge.width = 1,
    edge.color = "grey50",
    edge.curved = 0.3,
    layout = l
  )
  dev.off()
  
}

