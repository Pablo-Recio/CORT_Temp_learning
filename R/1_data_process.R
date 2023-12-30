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
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
    mutate(trt = factor(trt,
     levels = c("A_23", "A_28", "B_23", "B_28"),
     labels=c("A_23"="Control Cold",
             "A_28"="Control Hot",
             "B_23"="CORT Cold",
             "B_28"="CORT Hot"))) %>%
data.frame()

data_reversal <- data %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(FC_reversal)) <= 15) %>%
   ungroup()  %>%
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
    mutate(trt = factor(trt,
     levels = c("A_23", "A_28", "B_23", "B_28"),
     labels=c("A_23"="Control Cold",
             "A_28"="Control Hot",
             "B_23"="CORT Cold",
             "B_28"="CORT Hot"))) %>%
    mutate(trial_reversal=as.numeric(trial_reversal)) %>%
  data.frame()

# Standarize data by trial (i.e. make the first trial where each individual participated their trial 1)
data_associative <- data_associative %>%
  group_by(lizard_id) %>%
  mutate(
    first_non_na = min(which(!is.na(FC_associative))),
    Associative_Trial = ifelse(!is.na(first_non_na),trial_associative - first_non_na + 1, trial_associative))%>%
    filter(Associative_Trial >= 1) %>%
  ungroup() %>% 
  data.frame()
  
write.csv(data_associative, file= "./output/databases_clean/data_associative.csv")
write.csv(data_reversal, file= "./output/databases_clean/data_reversal.csv")
