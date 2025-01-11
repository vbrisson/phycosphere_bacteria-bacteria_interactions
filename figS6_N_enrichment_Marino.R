## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")
library(FSA)
source('utils.R')


permil_df <- read.csv("data/SIP_permil_v2.csv", check.names = FALSE)
sample_info <- read.csv("data/SIP_sample_info.csv")
sample_info$microplate <- as.character(sample_info$microplate)

permil_append <- append_xnet(permil_df, sample_info)

permil_app_stat <- permil_append %>%
  group_by(treatment, ring, strain, microplate) %>%
  summarize(N_q25 = quantile(N_permil, probs = 0.25),
            N_q50 = quantile(N_permil, probs = 0.5),
            N_q75 = quantile(N_permil, probs = 0.75),
            n = n()
            # cnet_mean = mean(Cnet),
            # cnet_sd = sd(Cnet),
            # replaced because we look for total incorp rate
  )


# global analysis
ggplot() +
  geom_sina(data = permil_append,
            aes(x=treatment, y=N_permil, color=treatment),
            maxwidth = 0.5,
            alpha=0.3,
            size=0.6) +
  geom_errorbar(data = permil_app_stat,
                aes(x=treatment, ymin=N_q50, ymax=N_q50),
                width = 0.3, size=0.8, color='black') +
  geom_errorbar(data = permil_app_stat,
                aes(x=treatment, ymin=N_q25, ymax=N_q75),
                width = 0.15, size=0.4, color='black') +
  geom_text(data = permil_app_stat,
            aes(x=treatment, y=-100, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  scale_y_continuous(breaks = append(seq(0, 800, 200), seq(1000, 5000, 1000))) +
  scale_y_break(c(900, 905), scales=0.4) +
  facet_grid(microplate ~ ring) +
  # scale_y_continuous(breaks = seq(0, 5000, 250)) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  # Alcani, Devosi, Marino, none
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        text = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = "black", size = 8),
        axis.ticks = element_line(colour = 'black', size=0.2),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.y = element_blank()
  )

ggsave("SIP_N_permil.pdf", width = 5, height = 6)


# # statistical test
# df_outer <- df_vis[df_vis$distance=='outer',]
# df_outer['treatment_p'] <- df_outer['treatment']!='none'
# kruskal.test(value~treatment_p, data=df_outer)
# 
# kruskal.test(value~treatment, data=df_outer)
# dunnTest(value~treatment, data=df_outer, method='holm')
# 
# t.test(df_vis[df_vis$treatment=='Alcanivorax',4], df_vis[df_vis$treatment=='Devosia',4])
# t.test(df_vis[df_vis$treatment=='Alcanivorax',4], df_vis[df_vis$treatment=='Marinobacter',4])
# t.test(df_vis[df_vis$treatment=='Devosia',4], df_vis[df_vis$treatment=='Marinobacter',4])
# 
# t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Alcanivorax',4])
# t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Devosia',4])
# t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Marinobacter',4])