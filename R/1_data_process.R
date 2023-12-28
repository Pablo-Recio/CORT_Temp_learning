####################################
# 1_data_process
####################################

# Packages
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, forcats)
source(here("R", "func.R"))

# Load data
data  <-  read.csv("./data/Learning.csv")

# Remove individuals who did not participate (more than 15 NAs), remove trials [36-40] (Only in associative) and split treatment into Temp and Cort
data_associative <- data %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(FC_associative)) <= 15) %>%
    filter(trial_associative <= 35) %>%
  ungroup()  %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
    mutate(temp = factor(temp,
     levels = c("23", "28"),
     labels=c("23"="Cold", "28"="Hot"))) %>%
    mutate(cort = factor(cort,
     levels = c("A", "B"),
     labels=c("A"="Control", "B"="CORT"))) %>%
  data.frame()

data_reversal <- data %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(FC_reversal)) <= 15) %>%
   ungroup()  %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
    mutate(temp = factor(temp,
     levels = c("23", "28"),
     labels=c("23"="Cold", "28"="Hot"))) %>%
    mutate(cort = factor(cort,
     levels = c("A", "B"),
     labels=c("A"="Control", "B"="CORT"))) %>%
    mutate(trial_reversal = ordered(trial_reversal)) %>%
  data.frame()

# Standarize data by trial (i.e. make the first trial where each individual participated their trial 1)
data_associative <- data_associative %>%
  group_by(lizard_id) %>%
  mutate(
    first_non_na = min(which(!is.na(FC_associative))),
    Associative_Trial = ifelse(!is.na(first_non_na),trial_associative - first_non_na + 1, trial_associative))%>%
    filter(Associative_Trial >= 1) %>%
    mutate(Associative_Trial = fct_inorder(as.factor(Associative_Trial))) %>%
  ungroup() %>% 
  data.frame()
  
write.csv(data_associative, file= "./output/Checking/data_associative.csv")

# Split data by species
deli_associative <- data_associative %>% filter(species == "delicata")
write.csv(deli_associative, file= "./output/databases_clean/deli_associative.csv")
guich_associative <- data_associative %>% filter(species == "guichenoti")
write.csv(guich_associative, file= "./output/databases_clean/guich_associative.csv")

# Split data by species
deli_reversal <- data_reversal %>% filter(species == "delicata")
write.csv(deli_reversal, file= "./output/databases_clean/deli_reversal.csv")
guich_reversal <- data_reversal %>% filter(species == "guichenoti")
write.csv(guich_reversal, file= "./output/databases_clean/guich_reversal.csv")
