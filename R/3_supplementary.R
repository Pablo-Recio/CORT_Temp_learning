####################################
# 3_Supplementary material
####################################

pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom)
source(here("R", "func.R"))

## Histogram ages by treatment
hist<-ggplot(data, aes(x = age.start, y = after_stat(density))) +  
  geom_density(alpha = 0.7,
  aes(fill = interaction(cort, temp))) +
  scale_fill_manual(values = c("CORT.Cold"="darkblue", "Control.Cold"="cyan", "CORT.Hot"="black", "Control.Hot"="grey"),
    labels=c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
  ) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +  
  theme(strip.background = element_blank()) +
  labs(y="Density", x = "Age", fill="Treatments") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),  
    axis.text = element_text(size = 10, family = "sans"),  
    legend.title = element_text(size = 14, family = "sans"), 
    legend.text = element_text(size = 12, family = "sans")  
  )
print(hist)
ggsave("./output/figures/hist.png", plot=hist, width = 30, height = 12, units = "cm", dpi = 1000) 