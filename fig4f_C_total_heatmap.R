library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")
library(FSA)
library(scales)
library(tidyverse)
source("utils.R")


## Combines count, Cnet, biovolume (day 14) for fig 4f
get_total_incorp_df <- function(){
  # import data
  count_bact <- read.csv("data/flow2_bact.csv")
  count_bact <- subset(count_bact, Time==14)
  count_bact$Ring[count_bact$Ring==1] <- 'inner'
  count_bact$Ring[count_bact$Ring==2] <- 'outer'
  count_pt <- read.csv("data/flow2_pt.csv")
  count_pt <- subset(count_pt, Time==14)
  cnet <- read.csv("data/SIP_cnet_v2.csv")
  cnet$Cnet <- cnet$Cnet*100 # to percent
  cnet_info <- read.csv("data/SIP_sample_info.csv")
  mass <- read.table("data/mass.csv", header = TRUE, check.names = F, sep = ',', nrows = 3)  # [fg C]
  doc_microplate <- read.csv('data/DesignResults/DOC_pm_alph_0.34_N_6.csv', header = TRUE)
  
  count_pt$Treatment[count_pt$Treatment==1] <- 'Devosia'
  count_pt$Treatment[count_pt$Treatment==2] <- 'none'
  count_pt$Treatment[count_pt$Treatment==3] <- 'Alcanivorax'
  count_pt$Treatment[count_pt$Treatment==4] <- 'Marinobacter'
  
  # Here we assume that bacterial cell number and Cnet are conditionally 
  # independent given algal cell number. For each microplate, the amount of 
  # bacterial incorporation of algal carbon is estimated by multiplying Cnet 
  # average with total bacterial number per ring. Optionally, it can then be 
  # divided by algal cell number to disentangle the algal effect. Normalized 
  # values can be used for statistical analysis, e.g., average or sd.
  
  # summarize count data
  count_bact_summ <- count_bact %>%
    group_by(Treatment, Ring) %>%
    summarize(abd_mean = mean(Abundance), 
              abd_sd = sd(Abundance),
              abd_med = median(Abundance),
              count_total = mean(Abundance) * 6, # per microplate
              n_count = n())
  count_bact_summ$Treatment[count_bact_summ$Treatment==1] <- 'Devosia'
  count_bact_summ$Treatment[count_bact_summ$Treatment==2] <- 'none'
  count_bact_summ$Treatment[count_bact_summ$Treatment==3] <- 'Alcanivorax'
  count_bact_summ$Treatment[count_bact_summ$Treatment==4] <- 'Marinobacter'
  
  count_pt_summ <- count_pt %>%
    group_by(Treatment) %>%
    summarize(Abundance_mean = mean(Abundance),
              Abundance_sd = sd(Abundance),
              n = n())
  
  # summarize cnet
  cnet$microplate <- NA
  cnet$ring <- NA
  cnet$treatment <- NA
  cnet$strain <- NA
  for (s in unique(cnet_info$sample_name)){
    m <- cnet_info$microplate[cnet_info$sample_name==s]
    r <- cnet_info$ring[cnet_info$sample_name==s]
    t <- cnet_info$treatment[cnet_info$sample_name==s]
    st <- cnet_info$strain[cnet_info$sample_name==s]
    
    cnet$microplate[cnet$sample_name==s] <- m
    cnet$ring[cnet$sample_name==s] <- r
    cnet$treatment[cnet$sample_name==s] <- t
    cnet$strain[cnet$sample_name==s] <- st
  }
  
  cnet_summ <- cnet %>%
    group_by(treatment, ring, strain) %>%
    summarize(cnet_mean = mean(Cnet),
              cnet_q50 = median(Cnet),
              cnet_sd = sd(Cnet),
              n_cnet = n()
    )
  
  # merge cell count, Cnet, single cell carbon mass, microplate DOC estimation
  # and algal cell count
  df <- select(cnet_summ, treatment, ring, strain, cnet_mean, cnet_q50, n_cnet)
  df$abd_mean <- NA
  df$n_count <- NA
  df$mass_c <- NA
  df$spatial <- NA
  df$count_pt <- NA
  df <- df[df$strain!='none',]
  
  for (row in 1:nrow(df)){
    t <- df$treatment[row]
    r <- df$ring[row]
    st <- df$strain[row]
    df[row, 'abd_mean'] <- count_bact_summ$abd_mean[
      count_bact_summ$Treatment == t & count_bact_summ$Ring == r
    ]
    df[row, 'n_count'] <- count_bact_summ$n_count[
      count_bact_summ$Treatment == t & count_bact_summ$Ring == r
    ]
    df[row, 'mass_c'] <- mass$C_mass_fg[mass$Genus==st]
    df[row, 'spatial'] <- ifelse(r=='inner', 
                                 doc_microplate[14,2] / doc_microplate[14,1], 
                                 doc_microplate[14,3] / doc_microplate[14,1])
    df[row, 'count_pt'] <- count_pt_summ$Abundance_mean[count_pt_summ$Treatment==t]
  }
  
  ## Display by estimate of total carbon
  # culture vol 70*6 ul, corr factor 0.75, total incubation 14 d
  # corr factor accounts for the volume halving on day 5
  df$c_incorp <- df$cnet_mean/100 * df$abd_mean * 0.07*6 * 0.75 * df$mass_c * 1e-6
  # [%] [cells ml-1] [ml] [] [fg C cell-1] [1e-6 ng fg-1] = [ng C]
  
  # Merge c_incorp across inner/outer wells
  df$c_incorp_total <- 0
  for (t in unique(df$treatment)){
    df$c_incorp_total[df$treatment==t] <- sum(df$c_incorp[df$treatment==t])
  }
  
  return(df)
}


## Main
cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet <- cnet[cnet$Cnet>0,]
cnet$Cnet <- cnet$Cnet*100 # to percent
cnet_info <- read.csv("data/SIP_sample_info.csv")

df <- get_total_incorp_df()
df$trt_loc <- paste(df$treatment, df$ring, sep='_')
df <- df[, c('trt_loc', 'c_incorp_total', 'c_incorp', 'spatial', 'cnet_q50', 
             'abd_mean', 'mass_c')]

df_vis <- 
  sweep(as.matrix(df[,2:7]), 2, as.matrix(df[7,2:7]), '/') %>% 
  as_tibble() %>%
  rownames_to_column("trt_loc") %>%
  pivot_longer(-trt_loc, names_to = "v2", values_to = "value")

df_long <- 
  as.matrix(df[,2:7]) %>% 
  as_tibble() %>%
  rownames_to_column("trt_loc") %>%
  pivot_longer(-trt_loc, names_to = "v2", values_to = "value")

df_vis$trt_loc <- factor(df_vis$trt_loc, levels=c(7,1:6))
df_vis$v2 <- factor(df_vis$v2, levels=rev(c('c_incorp_total', 'c_incorp', 'spatial',
                                            'cnet_q50', 'abd_mean', 'mass_c')))

## Heatmap Plot
ggplot() +
  geom_tile(data = df_vis, aes(trt_loc, v2, fill = value), color = "transparent") +
  geom_text(data = df_long, aes(trt_loc, v2, label = sprintf("%0.2f", round(value, digits = 2))), color = "black", size = 2) +
  # scale_fill_gradient( trans = 'log2' ) +
  scale_fill_distiller(palette = "PRGn", trans='log2', direction=1) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size=0.2),
        strip.background = element_blank(),
        strip.text = element_blank()
  )

ggsave("C_incorp_total_v3.pdf", width = 2.8, height = 1.9)

## plot, average over microplates
df_summ %>% 
  # arrange(c_incorp_mean) %>%
  # mutate(SDPos = cumsum(c_incorp_mean)) %>%
  mutate(treatment = factor(treatment, levels = c('none', 'Alcanivorax','Marinobacter', 'Devosia'))) %>%
  ggplot(aes(x = treatment, y = c_incorp_mean, 
             fill = treatment, alpha=factor(ring, levels=c("outer","inner")))) +
  geom_bar(stat = "identity", color='black', size=0.2, width=0.7) +
  geom_errorbar(aes(ymin = c_incorp_mean-c_incorp_sd, ymax = c_incorp_mean+c_incorp_sd), 
                width=0.1, linewidth=0.2, position = "identity", color='black',
                alpha =1) +
  scale_y_continuous(breaks= seq(0,80,5)) +
  scale_alpha_discrete(range=c(0.3, 1)) +
  scale_fill_manual(values=c("white", "#bd2529", "#1b6b3b", "#4159a7")) +  
  coord_flip() + 
  facet_wrap(~ring, ncol=2) +
  ylab('Carbon mass incorporation (ng C)') +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(colour = 'black', size=0.2),
        strip.background = element_blank(),
        strip.text = element_blank()
  )

ggsave("C_incorp_total_v2.pdf", width = 3.0, height = 1.7)

write.csv(df_summ, "data/c_mass_incorp.csv")