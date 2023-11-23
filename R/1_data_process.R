####################################
# 1_data_process
####################################

## Load packages
	install.packages("pacman")
	pacman::p_load(tidyverse, readxl, ggplot2, dplyr, magrittr)

## Load data	
data<-read.csv("./data/learning.csv")

data<- data %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>% 
    mutate(temp = factor(temp,
     levels = c("23", "28"), 
     labels=c("23"="Cold", "28"="Hot"))) %>%
    mutate(cort = factor(cort,
     levels = c("A", "B"), 
     labels=c("A"="Control", "B"="CORT"))) %>%
	mutate(species = factor(species,
	 labels=c("delicata"="L. delicata", "guichenoti"="L. guichenoti")))
  data.frame()


## Histogram ages by treatment
hist<-ggplot(data, aes(x = age.start, y = ..density..)) +  
  geom_density(alpha = 0.7,
  aes(fill = interaction(cort, temp))) +
  scale_fill_manual(values = c("CORT.Cold"="darkblue", "Control.Cold"="cyan", "CORT.Hot"="black", "Control.Hot"="grey"),
    labels=c("CORT-Cold (n=11)", "Control-Cold (n=12)", "CORT-Hot (n=11)", "Control-Hot (n=10)")
  ) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +  
  theme(strip.background = element_blank()) +
  labs(y="", x = "Age") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),  
    axis.text = element_text(size = 10, family = "sans"),  
    legend.title = element_text(size = 14, family = "sans"), 
    legend.text = element_text(size = 12, family = "sans")  
  )
print(hist)
ggsave("./output/figures/hist.png", plot=hist, width = 30, height = 12, units = "cm", dpi = 1000) 

## Merge



# At the end of processing you write the file
	write_csv(data_activity, "./output/data/data_activity.csv")

