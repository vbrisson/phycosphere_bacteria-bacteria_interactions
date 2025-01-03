library('ggplot2')
library('ggforce')
library("dplyr")
source("utils.R")  # Continued from "fig4_pm_count.R"
# detach(package:Rmisc)
# detach(package:plyr)


df <- get_df()

df_stat <- df %>%
  group_by(Treatment, Ring, Time) %>%
  summarize(q25 = quantile(Abundance, probs = 0.25), 
            q50 = quantile(Abundance, probs = 0.5),
            q75 = quantile(Abundance, probs = 0.75),
            n = n()
            )


##### PLOT ABUNDANCE, DAY 14
## OUTER RING
df_vis <- df[df$Ring==2 & df$Time==14,]
df_stat_vis <- df_stat[df_stat$Ring==2 & df_stat$Time==14,]

ggplot() +
  geom_sina(data = df_vis, 
            aes(Treatment, Abundance, color=Treatment),
            scale = 'width',
            size=1.2,
            maxwidth = 0.45, shape=21, fill="white", stroke=0.5) +
  geom_errorbar(data=df_stat_vis, 
                aes(x=Treatment, ymin=q25, ymax=q75),
                width = 0.15, size=0.3, color='black') + 
  geom_errorbar(data=df_stat_vis, 
                aes(x=Treatment, ymin=q50, ymax=q50),
                width = 0.3, size=0.6, color='black') + 
  scale_y_continuous(breaks= seq(0,1.4e7,4e6)) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
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

ggsave("fig4_abundance_outer.pdf", width = 1.5, height = 2.5)


## stats test outer ring
kruskal.test(Abundance ~ Treatment, data = df_vis)
pairwise.wilcox.test(
  df_vis$Abundance, 
  df_vis$Treatment, 
  p.adjust.method = "BH"
)