library("ggplot2")
library(ggforce)
library("dplyr")
source('utils.R')

permil_df <- read.csv("data/SIP_permil_v2.csv", check.names = FALSE)
sample_info <- read.csv("data/SIP_sample_info.csv")
sample_info$microplate <- as.character(sample_info$microplate)

## add sample details to df using sample_info
permil_append <- append_xnet(permil_df, sample_info)


## summarize df
permil_df_stat <- permil_append %>%
  group_by(ring, treatment, strain) %>%
  summarize(q25 = quantile(ROIAREA, probs = 0.25), 
            q50 = quantile(ROIAREA, probs = 0.5),
            q75 = quantile(ROIAREA, probs = 0.75),
            mean = mean(ROIAREA),
            n = n()
            )


## plot outer strains
permil_append_vis <- subset(permil_append, ring=='outer')
permil_df_stat_vis <- subset(permil_df_stat, ring=='outer')

ggplot() +
  geom_sina(data = permil_append_vis,
            aes(x=treatment, y=ROIAREA, color=treatment),
            maxwidth = 0.5,
            alpha=0.3,
            size=0.6) +
  geom_errorbar(data = permil_df_stat_vis,
                aes(x=treatment, ymin=q50, ymax=q50),
                width = 0.3, size=0.8, color='black') +
  geom_errorbar(data = permil_df_stat_vis,
                aes(x=treatment, ymin=q25, ymax=q75),
                width = 0.15, size=0.4, color='black') +
  geom_text(data = permil_df_stat_vis,
            aes(x=treatment, y=0, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  scale_y_continuous(limits = c(0,4)) +
  ylab('Cell area (um2)') +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 8),
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.ticks = element_line(colour = 'black', size=0.2),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        text = element_text(size = 8)
  )

ggsave("figures/figS_size_outer.pdf", width = 3, height = 3)


## plot inner strains
permil_append_vis <- subset(permil_append, ring=='inner')
permil_df_stat_vis <- subset(permil_df_stat, ring=='inner')

ggplot() +
  geom_sina(data = permil_append_vis,
            aes(x=treatment, y=ROIAREA, color=treatment),
            maxwidth = 0.5,
            alpha=0.3,
            size=0.6) +
  geom_errorbar(data = permil_df_stat_vis,
                aes(x=treatment, ymin=q50, ymax=q50),
                width = 0.3, size=0.8, color='black') +
  geom_errorbar(data = permil_df_stat_vis,
                aes(x=treatment, ymin=q25, ymax=q75),
                width = 0.15, size=0.4, color='black') +
  geom_text(data = permil_df_stat_vis,
            aes(x=treatment, y=-0.2, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  # scale_y_continuous(limits = c(0,5)) +
  ylab('Cell area (um2)') +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 8),
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.ticks = element_line(colour = 'black', size=0.2),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        text = element_text(size = 8)
  )

ggsave("figures/figS_size_inner.pdf", width = 2.5, height = 3)



## statistical test
permil_append_vis <- subset(permil_append, ring=='outer')
permil_df_stat_vis <- subset(permil_df_stat, ring=='outer')

kruskal.test(ROIAREA ~ treatment, data = permil_append_vis)
pairwise.wilcox.test(permil_append_vis$ROIAREA, permil_append_vis$treatment, p.adjust.method = "BH")