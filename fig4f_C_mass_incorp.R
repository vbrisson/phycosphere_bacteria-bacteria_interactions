source("utils.R")

df_summ <- get_total_incorp_df()$df_summ

# plot, average over microplates

df_summ %>% mutate(SD = c_incorp_sd) %>%
  group_by(treatment) %>%
  mutate(SDPos = cumsum(c_incorp_mean)) %>%
  ggplot(aes(x = reorder(treatment, -c_incorp_mean), y = c_incorp_mean, fill = treatment,
             alpha=factor(ring, levels=c("outer","inner"))
             )) +
  geom_bar(stat = "identity", color='black', size=0.2, width=0.7) +
  # geom_text(aes(label = paste(round(c_incorp_mean,3), '±', round(c_incorp_sd,3))), 
  #           position = position_stack(vjust =  0.5)) +
  geom_errorbar(aes(ymin = SDPos-SD, ymax = SDPos+SD), 
                width=0.1, 
                linewidth=0.2, 
                position = "identity", 
                color='black',
                alpha =1) +
  scale_alpha_discrete(range=c(0.3, 1)) +
  # scale_fill_manual(values=c("#FFDFE1", "#D2DCFB", "#F8DFC4", "#D1D3D4")) +
  scale_fill_manual(values=c("#bd2529", "#4159a7", "#1b6b3b", "white")) +  
  # Alcani, Devosi, Marino, none
  ylab('carbon mass incorporation (ng C)') +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(colour = 'black', size=0.2),
  )

ggsave("figures/C_incorp_total.pdf", width = 1.6, height = 1.7)