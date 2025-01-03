source('utils.R')
# detach(package:Rmisc)
# detach(package:plyr)
library("dplyr")

df <- get_df()
df_alg <- get_df_alg()
df_fold <- get_df_fold(df, df_alg)
df <- df[df$Strain!='None',]
df$Ring[df$Ring==1] <- 'inner'
df$Ring[df$Ring==2] <- 'outer (Marinobacter)'


# bacteria abundance stat
df_stat <- df %>%
  group_by(Treatment, Ring, Time) %>%
  summarize(q25 = quantile(Abundance, probs = 0.25), 
            q50 = quantile(Abundance, probs = 0.5),
            q75 = quantile(Abundance, probs = 0.75),
            n = n()
            )


### Bacterial Abundance by Time
ggplot() +
  geom_point(data = df, 
            aes(x=Time, y=Abundance/1e6, color=Treatment), size=0.8) +
  facet_grid(cols = vars(Treatment), rows = vars(Ring)) +
  geom_line(data=df_stat, aes(x=Time, y=q50/1e6, color=Treatment), size=0.4) +
  geom_errorbar(data=df_stat, 
                aes(x=Time, ymin=q25/1e6, ymax=q75/1e6),
                width = 0.5, size=0.2, color='black') + 
  geom_errorbar(data=df_stat, 
                aes(x=Time, ymin=q50/1e6, ymax=q50/1e6),
                width = 1, size=0.3, color='black') + 
  geom_text(data=df_stat,
            aes(x=Time, y=0, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  labs(y = "Bacteria (million ml-1)", x = "Time (d)") +
  theme(strip.background = element_blank(),
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

ggsave("figS_bact_abundance_global.pdf", width = 6.5, height = 4)