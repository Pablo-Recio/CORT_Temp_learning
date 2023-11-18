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

pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom)

source(here("R", "func.R"))
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
 
# Load data
	data  <-  read.csv("./data/Learning.csv")
#
#
#
#| label: dataclean
#| echo: false
#| warning: false

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
  data.frame() 


# Standarize data by trial (i.e. make the first trial where each individual participated their trial 1)
data_associative <- data_associative %>%
  group_by(lizard_id) %>%
  mutate(
    first_non_na = min(which(!is.na(FC_associative))),
    Associative_Trial = ifelse(!is.na(first_non_na),trial_associative - first_non_na + 1, trial_associative))%>%
    filter(Associative_Trial >= 1) %>%
   ungroup() %>% data.frame()
  
write.csv(data_associative, file= "./output/Checking/data_associative.csv")

# Split data by species
deli_associative <- data_associative %>% filter (species == "delicata")
guich_associative <- data_associative %>% filter (species == "guichenoti")

# Split data by species
deli_reversal <- data_reversal %>% filter (species == "delicata")
guich_reversal <- data_reversal %>% filter (species == "guichenoti")
#
#
#
#
#| label: models1
#| echo: false
#| warning: false

# Generalized Mixed Models delicata
GMM_deli <- glmer(FC_associative ~ Associative_Trial*temp*cort + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = deli_associative, family = binomial(link = "logit"))
summary(GMM_deli)

emmeans(GMM_deli, pairwise ~ temp*Associative_Trial)
emmeans(GMM_deli, pairwise ~ cort*Associative_Trial)

# Generalized Mixed Models guichenoti
GMM_guich <- glmer(FC_associative ~ Associative_Trial*temp*cort + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = guich_associative, family = binomial(link = "logit"))
summary(GMM_guich)

emmeans(GMM_guich, pairwise ~ temp*Associative_Trial)
emmeans(GMM_guich, pairwise ~ cort*Associative_Trial)
emmeans(GMM_guich, pairwise ~ temp*cort*Associative_Trial)

#
#
#
#| label: modelsasso
#| echo: false
#| warning: false

# Bayesian model delicata
  deli_mod <- brm(FC_associative ~ cort*temp*Associative_Trial + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = deli_associative, family = bernoulli(link = "logit"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(deli_mod)
  plot(deli_mod)

# Bayesian model guichenoti
  guich_mod <- brm(FC_associative ~ cort*temp*Associative_Trial + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = guich_associative, family = bernoulli(link = "logit"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(guich_mod)
  plot(guich_mod)
# R2
bayes_R2(deli_mod)
bayes_R2(guich_mod)

# Extract posteriors
posterior <- posterior_samples(mod, pars = "^b")
plyr::ldply(lapply(posterior, mean))
plyr::ldply(lapply(posterior, function(x) quantile(x, c(0.025, 0.975))))
temp23 <- c(posterior[,"b_Intercept"], (posterior[,"b_Intercept"]+ posterior[,"b_trtB_23"]))
temp28 <- c((posterior[,"b_Intercept"] + posterior[,"b_trtB_28"]), posterior[,"b_trtB_23"])

pmcmc(temp23 - temp28)
dharma(GMM_deli2)


#
#
#
#| label: modelsrev
#| echo: false
#| warning: false

# Bayesian model delicata
  deli_rev <- brm(FC_reversal ~ cort*temp*trial_reversal + group*trial_reversal + (1 + trial_reversal|lizard_id), data = deli_reversal, family = bernoulli(link = "logit"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(deli_rev)
  plot(deli_rev)

# Bayesian model guichenoti
  guich_rev <- brm(FC_reversal ~ cort*temp*trial_reversal + group*trial_reversal + (1 + trial_reversal|lizard_id), data = guich_reversal, family = bernoulli(link = "logit"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(guich_rev)
  plot(guich_rev)
# R2
bayes_R2(deli_rev)
bayes_R2(guich_rev)

# Extract posteriors
posterior <- posterior_samples(mod, pars = "^b")
plyr::ldply(lapply(posterior, mean))
plyr::ldply(lapply(posterior, function(x) quantile(x, c(0.025, 0.975))))
temp23 <- c(posterior[,"b_Intercept"], (posterior[,"b_Intercept"]+ posterior[,"b_trtB_23"]))
temp28 <- c((posterior[,"b_Intercept"] + posterior[,"b_trtB_28"]), posterior[,"b_trtB_23"])

pmcmc(temp23 - temp28)
dharma(GMM_deli2)


#
#
#
#| label: plots
#| echo: false
#| warning: false
 
# Plot delicata associative
plot_deli <- ggplot(deli_associative, aes(x = Associative_Trial, y = FC_associative, color=interaction(cort,temp))) +  
  geom_smooth(method = "glm", formula = y ~ x, method.args = list(family = "binomial")) +
  facet_grid(as.factor(deli_associative$group)) + 
  labs(y = ""Probability of correct choice"", x = "Trial")  
print(plot_deli)

# Plot guichenoti associative
plot_guich <- ggplot(guich_associative, aes(x = Associative_Trial, y = FC_associative, color=interaction(cort,temp))) +  
  geom_smooth(method = "glm", formula = y ~ x, method.args = list(family = "binomial"),na.action = na.exclude) +
  facet_grid(as.factor(guich_associative$group)) + 
  labs(y = "Probability of correct choice", x = "Trial") 
print(plot_guich)


# Plot delicata reversal

# Plot guichenoti reversal

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
