#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: setup

pacman::p_load(tidyverse, flextable)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: dataload
#| echo: false
#| warning: false
 
# Load packages
pacman::p_load(tidyverse, flextable, brms, ggplot2)

# Load data
	data  <-  read.csv("./data/Learning.csv")
#
#
#
#| label: dataclean
#| echo: false
#| warning: false

# Removing individuals who did not participate (more than 15 NAs)
data_associative <- data %>%
  group_by(lizard_id) %>%
  filter(sum(is.na(FC_associative)) <= 15) %>%
  ungroup()
data_associative <- as.data.frame(data_associative)

# Removing trials [36-40] (Only in associative)
data_associative <- data_associative %>%
group_by(lizard_id) %>%
filter(trial_associative <= 35) %>%
ungroup()
data_associative <- as.data.frame(data_associative)

# Standarize data by trial (i.e. make the first trial where each individual participated their trial 1)

data_associative <- data_associative %>%
  group_by(lizard_id) %>%
  mutate(
    first_non_na = min(which(!is.na(FC_associative))),
    Associative_Trial = ifelse(!is.na(first_non_na),trial_associative - first_non_na + 1, trial_associative))%>%
    filter(Associative_Trial >= 1) %>%
   ungroup()
  
data_associative <- as.data.frame(data_associative)
write.csv(data_associative, file= "./output/Checking/data_associative.csv")

# Split data
deli_associative <- data_associative %>% filter (species == "delicata")
guich_associative <- data_associative %>% filter (species == "guichenoti")

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: fig-fig1
#| fig-cap: "Activity budget of the focal group of chacma baboons (Papio ursinus) at the Cape Peninsula, South Africa, during the study period. The data are presented as the mean percentage of time spent in each activity category per hour (± SE)."

data %>% ggplot(aes(x = Lizard_id))  + geom_histogram()


#
#
#
#
#
#| label: tbl-tb1
#| tbl-cap: "Summary of the data"

tab <- data %>% 
  group_by(Lizard_id) %>% 
  summarise(
    n = n()
  ) 
flextable(tab[2:3,])
#
#
#
#
#
