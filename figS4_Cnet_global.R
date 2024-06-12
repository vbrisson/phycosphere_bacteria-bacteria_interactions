## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")
library(FSA)
source('utils.R')


cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet <- cnet[cnet$Cnet>0,]
cnet$Cnet <- cnet$Cnet*100 # to percent
cnet_info <- read.csv("data/SIP_sample_info.csv")
cnet_append <- append_xnet(cnet, cnet_info)
cnet_append_stat <- cnet_append %>%
  group_by(treatment, ring, microplate) %>%
  summarize(q25 = quantile(Cnet, probs = 0.25),
            q50 = quantile(Cnet, probs = 0.5),
            q75 = quantile(Cnet, probs = 0.75),
            n = n(),
            max = max(Cnet)
  )


# All conditions
ggplot() +
  geom_sina(data = cnet_append,
            aes(x=ring, y=Cnet, color=treatment),
            maxwidth = 0.5,
            alpha=0.3,
            size=0.6) +
  facet_grid(microplate ~ treatment) +
  geom_errorbar(data = cnet_append_stat,
                aes(x=ring, ymin=q50, ymax=q50),
                width = 0.3, size=0.8, color='black') +
  geom_errorbar(data = cnet_append_stat,
                aes(x=ring, ymin=q25, ymax=q75),
                width = 0.15, size=0.4, color='black') +
  geom_text(data = cnet_append_stat,
            aes(x=ring, y=-0.2, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  # scale_y_continuous(breaks = append(seq(0, 9, 3), seq(10, 25, 10))) +
  # scale_y_break(c(9, 10), scales=0.1) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  # Alcani, Devosi, Marino, none
  theme(strip.background = element_rect(fill=NA),
         panel.background = element_rect(fill = "transparent", color = NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.background = element_rect(fill = "transparent", color = NA),
         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
         legend.position = "none",
         text = element_text(colour = 'black', size = 8),
         axis.text = element_text(colour = "black", size = 8),
         axis.ticks = element_line(colour = 'black', size=0.2),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank(),
         axis.line.y.right = element_blank()
  )

ggsave("figures/SIP_Cnet_global.pdf", width = 6, height = 5)