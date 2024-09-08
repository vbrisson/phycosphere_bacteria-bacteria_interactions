## Combines count and Cnet data in one plot (day 14)
library("dplyr")
library(ggplot2)
library(scales)
source('utils.R')


# import and clean data
count_bact <- read.csv("data/flow2_bact.csv")
count_bact <- subset(count_bact, Time==14)
count_bact$Ring[count_bact$Ring==1] <- 'inner'
count_bact$Ring[count_bact$Ring==2] <- 'outer'
count_bact$Treatment[count_bact$Treatment==1] <- 'Devosia'
count_bact$Treatment[count_bact$Treatment==2] <- 'none'
count_bact$Treatment[count_bact$Treatment==3] <- 'Alcanivorax'
count_bact$Treatment[count_bact$Treatment==4] <- 'Marinobacter'
count_pt <- read.csv("data/flow2_pt.csv")
count_pt <- subset(count_pt, Time==14)
count_pt$Treatment[count_pt$Treatment==1] <- 'Devosia'
count_pt$Treatment[count_pt$Treatment==2] <- 'none'
count_pt$Treatment[count_pt$Treatment==3] <- 'Alcanivorax'
count_pt$Treatment[count_pt$Treatment==4] <- 'Marinobacter'
cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet <- cnet[cnet$Cnet>0,]
cnet_info <- read.csv("data/SIP_sample_info.csv")

# summarize cnet
cnet_stat <- summarize_cnet(cnet, cnet_info)$df_summ

# get statistics from count data
count_bact_stat <- count_bact %>%
  group_by(Treatment, Microplate, Ring) %>%
  summarize(mean = mean(Abundance), 
            sd = sd(Abundance))

# merge count and cnet in to df
merged_stat <- merge_count_cnet(count_bact_stat, cnet_stat)


# plot inner ring
setEPS()
postscript("figures/fig5a_v2.eps", width = 1.6, height = 1.6)

ggplot(merged_stat[merged_stat$ring=='inner',], aes(x=cnet_q50, y=count_mean, colour=treatment)) + 
  geom_point(aes(size=n_cnet, fill=treatment, shape=ring), stroke=0.5) + 
  geom_errorbar(aes(ymax = count_mean+count_sd, ymin = count_mean-count_sd), width=0, linewidth=0.1) + 
  geom_errorbarh(aes(xmax = cnet_q75, xmin = cnet_q25), height=0, linewidth=0.1) + 
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans='log10', breaks=c(0.1)) +
  scale_y_continuous(trans='log10', breaks=c(10^6, 10^7)) +
  annotation_logticks(short = unit(0.1, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.2, "cm"),
                      size = 0.1) +
  scale_size(range = c(1, 2.5)) +
  scale_color_manual(values=c("#bd2529","#4159a7","#1b6b3b","#231f20")) +
  scale_fill_manual(values=c("#fcdcdf", "#d2daef", "#c4d2cc", "#ffffff")) +
    # Alcani, Devosi, Marino, none
  scale_shape_manual(values = c(16,21)) +  # Inner, outer
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_line(colour = "grey80", linewidth=0.2),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.15),
        axis.ticks = element_blank()
        )

dev.off()