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
  group_by(treatment, ring, strain) %>%
  summarize(N_q25 = quantile(N_permil, probs = 0.25),
            N_q50 = quantile(N_permil, probs = 0.5),
            N_q75 = quantile(N_permil, probs = 0.75),
            n = n()
            # cnet_mean = mean(Cnet),
            # cnet_sd = sd(Cnet),
            # replaced because we look for total incorp rate
  )


# Outer strain
permil_append_vis <- subset(permil_append, ring=='outer')
permil_app_stat_vis <- subset(permil_app_stat, ring=='outer')

ggplot() +
  geom_sina(data = permil_append_vis,
            aes(x=treatment, y=N_permil, color=treatment),
            maxwidth = 0.5,
            alpha=0.3,
            size=0.6) +
  geom_errorbar(data = permil_app_stat_vis,
                aes(x=treatment, ymin=N_q50, ymax=N_q50),
                width = 0.3, size=0.8, color='black') +
  geom_errorbar(data = permil_app_stat_vis,
                aes(x=treatment, ymin=N_q25, ymax=N_q75),
                width = 0.15, size=0.4, color='black') +
  geom_text(data = permil_app_stat_vis,
            aes(x=treatment, y=-0.2, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  scale_y_continuous(breaks = append(seq(0, 750, 250), seq(1000, 5000, 2000))) +
  scale_y_break(c(750, 751), scales=0.1) +
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
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.ticks = element_line(colour = 'black', size=0.2),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8)
  )

ggsave("figures/SIP_N_permil_outer.pdf", width = 1.5, height = 2.5)
# ggsave("figures/SIP_N_permil_outer_annot.pdf", width = 3, height = 3)


# statistical test
kruskal.test(N_permil ~ treatment, data = permil_append_vis)
pairwise.wilcox.test(
  permil_append_vis$N_permil, 
  permil_append_vis$treatment, 
  p.adjust.method = "BH"
)