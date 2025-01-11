library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")
library(FSA)
library(scales)
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
  cnet_info <- read.csv("data/SIP_sample_info.csv")
  mass <- read.table("data/mass.csv", header = TRUE, sep = ',', nrows=3)  # [fg]
  
  count_pt$Treatment[count_pt$Treatment==1] <- 'Devosia'
  count_pt$Treatment[count_pt$Treatment==2] <- 'none'
  count_pt$Treatment[count_pt$Treatment==3] <- 'Alcanivorax'
  count_pt$Treatment[count_pt$Treatment==4] <- 'Marinobacter'
  
  # Here we assume that bacterial cell number and Cnet are conditionally 
  # independent given algal cell number. For each microplate, the amount of 
  # bacterial incorporation of algal carbon is estimated by multiplying Cnet 
  # average with total bacterial number per ring. Normalized values can be used
  # for statistical analysis, e.g., average or sd.
  
  # summarize count data
  count_bact_summ <- count_bact %>%
    group_by(Treatment, Ring, Microplate) %>%
    summarize(abd_mean = mean(Abundance), 
              abd_sd = sd(Abundance),
              abd_med = median(Abundance),
              count_total = mean(Abundance) * 6
    ) # merge over 6 wells per ring
  count_bact_summ$Treatment[count_bact_summ$Treatment==1] <- 'Devosia'
  count_bact_summ$Treatment[count_bact_summ$Treatment==2] <- 'none'
  count_bact_summ$Treatment[count_bact_summ$Treatment==3] <- 'Alcanivorax'
  count_bact_summ$Treatment[count_bact_summ$Treatment==4] <- 'Marinobacter'
  
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
    group_by(treatment, microplate, ring, strain) %>%
    summarize(cnet_mean = mean(Cnet),
              cnet_sd = sd(Cnet),
              # replaced because we look for total incorp rate
              cnet_q50 = quantile(Cnet, probs = 0.5)
              # cnet_q25 = quantile(Cnet, probs = 0.25),
              # cnet_q75 = quantile(Cnet, probs = 0.75)
    )
  
  # merge count, cnet, mass, algal count in to df
  df <- select(cnet_summ, treatment, microplate, ring, strain, cnet_mean, cnet_q50)
  df$count_total <- NA
  df$mass_c <- NA
  df$count_pt <- NA
  df <- df[df$strain!='none',]
  
  for (row in 1:nrow(df)){
    t <- df$treatment[row]
    m <- df$microplate[row]
    r <- df$ring[row]
    st <- df$strain[row]
    df[row, 'count_total'] <- count_bact_summ$count_total[
      count_bact_summ$Treatment == t & 
        count_bact_summ$Ring == r &
        count_bact_summ$Microplate == m
    ]
    df[row, 'mass_c'] <- mass$CarbonMass_fg[mass$Genus==st]
    df[row, 'count_pt'] <- count_pt$Abundance[
      count_pt$Treatment==t & 
        count_pt$Microplate==m
    ]
  }
  
  ## Display by estimate of total carbon
  # culture vol 70 ul, corr factor 0.75 due to sampling at day 5, 
  # total incubation 14 d corr factor accounts for the volume halving on day 5
  df$c_incorp <- df$cnet_mean * df$count_total * 0.07 * 0.75 * df$mass_c * 1e-6
  # [%] [cells ml-1] [ml] [] [fg C cell-1] [1e-6 ng fg-1] = [ng C]
  
  df_summ <- df %>%
    group_by(treatment, ring) %>%
    summarize(
      # c_incorp_pt_mean = 100*mean(c_incorp_pt), 
      # c_incorp_pt_sd = 100*sd(c_incorp_pt),
      c_incorp_mean = mean(c_incorp),
      c_incorp_sd = sd(c_incorp)
    )
  
  return(df_summ)
}


cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet$Cnet <- cnet$Cnet*100 # to percent
cnet_info <- read.csv("data/SIP_sample_info.csv")

df_summ <- get_total_incorp_df()


## plot, average over microplates
df_summ[df_summ$treatment!='none',] %>% 
  mutate(SD = c_incorp_sd) %>%
  group_by(treatment) %>%
  mutate(SDPos = cumsum(c_incorp_mean)) %>%
  ggplot(aes(x = reorder(treatment, -c_incorp_mean), y = c_incorp_mean, fill = treatment,
             color = treatment, alpha=factor(ring, levels=c("outer","inner"))
             )) +
  geom_bar(stat = "identity", color='black', size=0.2, width=0.5) +
  # geom_text(aes(label = paste(round(c_incorp_mean,3), '±', round(c_incorp_sd,3))), 
  #           position = position_stack(vjust =  0.5)) +
  geom_errorbar(aes(ymin = SDPos-SD, ymax = SDPos+SD), 
                width=0.1, 
                linewidth=0.2, 
                position = "identity", 
                color='black',
                alpha =1) +
  scale_alpha_discrete(range=c(1, 0)) +
  # scale_fill_manual(values=c("#FFDFE1", "#D2DCFB", "#F8DFC4", "#D1D3D4")) +
  scale_color_manual(values=c("#bd2529", "#4159a7", "#1b6b3b", "white")) +  
  scale_fill_manual(values=c("#bd2529", "#4159a7", "#1b6b3b", "white")) +  
  # Alcani, Devosi, Marino, none
  ylab('carbon mass incorporation (ng C)') +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = NA),
        panel.border = element_blank(),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(colour = 'black', size=0.2),
  )

ggsave("C_incorp_total_v2.pdf", width = 1.8, height = 1.8)

write.csv(df_summ, "data/c_mass_incorp.csv")