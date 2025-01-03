## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
# detach(package:plyr)
library("ggbreak")
source('utils.R')


cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet <- cnet[cnet$Cnet>0,]
cnet$Cnet <- cnet$Cnet*100 # to percent
cnet_info <- read.csv("data/SIP_sample_info.csv")
cnet_append <- append_xnet(cnet, cnet_info)
cnet_outer <- cnet_append[(cnet_append$ring=='outer'),]
cnet_outer_stat <- cnet_outer %>%
   group_by(treatment, ring, microplate) %>%
   summarize(q25 = quantile(Cnet, probs = 0.25),
             q50 = quantile(Cnet, probs = 0.5),
             q75 = quantile(Cnet, probs = 0.75),
             n = n(),
             max = max(Cnet)
   )


# draw figures
ggplot() +
  geom_sina(data = cnet_outer, 
            aes(x=treatment, y=Cnet, color=treatment), 
            maxwidth = 0.8,
            alpha=0.4,
            size=0.3) +
  geom_errorbar(data=cnet_outer_stat, 
                aes(x=treatment, ymin=q25, ymax=q75),
                width = 0.15, color='black', size=0.3) + 
  geom_errorbar(data=cnet_outer_stat, 
                aes(x=treatment, ymin=q50, ymax=q50),
                width = 0.3, color='black', size=0.6) + 
  labs(x = "Influencer", y = "Cnet") +
  scale_y_continuous(breaks = append(seq(0, 9, 3), seq(10, 25, 10))) +
  scale_y_break(c(9, 10), scales=0.1) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +  
    # Alcani, Devosi, Marino, none
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 8),
        axis.ticks = element_line(colour = 'black', size=0.2),
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8)
        )

ggsave("figures/SIP_cnet_outer.pdf", width = 1.5, height = 2.5)
# ggsave("figures/SIP_cnet_day14_inner_break.pdf", width = 3, height = 4)


# statistical test
kruskal.test(Cnet ~ treatment, data = cnet_outer)
pairwise.wilcox.test(cnet_outer$Cnet, cnet_outer$treatment, p.adjust.method = "BH")

print(
  cat(
    'Alcanivorax:', (df_vis_stat$q50[1]-df_vis_stat$q50[4]) / df_vis_stat$q50[4], 
    '\n',
    'Devosia:', (df_vis_stat$q50[2]-df_vis_stat$q50[4]) / df_vis_stat$q50[4], 
    '\n',
    'Marinobacter:', (df_vis_stat$q50[3]-df_vis_stat$q50[4]) / df_vis_stat$q50[4],
    '\n'
    )
  )