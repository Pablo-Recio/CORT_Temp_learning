####################################
# 2_MODELS
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom)
#
######## 2.A) BUILD ALL MODELS AND EXTRACT POSTERIORS
#### Here we fit all the (brm) models and extract the posteriors per each of the treatments using output/clean_databases already processed using the code in "1_data_process.R"
#### All this posteriors are processed by fit_asso function and a glbal database with the extimates of Associative_Trial per each treatment level, species, and group is created
## Associative task
source(here("R", "func.R"))
species_levels <- unique(data_associative$species)
bias_levels <- unique(data_associative$group)
global_asso <- data.frame()
for(species_level in species_levels) {
  for(bias_level in bias_levels) {
    res <- fit_asso(sp = species_level, bias = bias_level)
    # Add columns for species and bias to the result
    res$species_level <- species_level
    res$bias_level <- bias_level
    # Combine the result with the global result
    global_asso <- bind_rows(global_asso, res)
  }
}
write.csv(global_asso, here("output/Checking/global_asso.csv"))
data_trial <- filter(data_associative, data_associative$species=="delicata")
data_trial$group <- relevel(data_trial$group, ref = "Blue")
data_trial$trt <- relevel(data_trial$trt, ref = "Control Cold")
m <- brm(FC_associative ~ trt*Associative_Trial + Associative_Trial*group + (1 + Associative_Trial|lizard_id),
                data = data_trial,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
posterior_trial <- as_draws(m)
write.csv(posterior_trial, here("output/Checking/posterior_trial.csv"))
data_trial_2 <- filter(data_associative, data_associative$species=="delicata")
data_trial_2$group <- relevel(data_trial$group, ref = "Blue")
data_trial_2$trt <- relevel(data_trial$trt, ref = "CORT Cold")
m2 <- brm(FC_associative ~ trt*Associative_Trial + Associative_Trial*group + (1 + Associative_Trial|lizard_id),
                data = data_trial_2,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
posterior_trial_2 <- as_draws(m2)
write.csv(posterior_trial_2, here("output/Checking/posterior_trial_2.csv"))
#
## Reversal task
source(here("R", "func.R"))
species_levels <- unique(data_reversal$species)
bias_levels <- unique(data_reversal$group)
global_rev <- data.frame()
for(species_level in species_levels) {
  for(bias_level in bias_levels) {
    res <- fit_rev(sp = species_levels, bias = bias_levels)
    # Add columns for species and bias to the result
    res$species_level <- species_level
    res$bias_level <- bias_level
    # Combine the result with the global result
    global_rev <- bind_rows(global_rev, res)
  }
}
write.csv(global_rev, here("output/Checking/global_rev.csv"))
#
######## 2.B) Make comparisons between estimates and create pmcmcs

# pmcmc for specific comparisons

# pmcmc for differences between 'Groups'

# pmcmc for differences between both species

######## 2.C) MAKE TABLE WITH POSTERIORS AND pmcmc

######## 2.E) PLOTS MODELS 
# Plot delicata associative
plot_deli <- ggplot(deli_associative, aes(x = Associative_Trial, y = FC_associative, color=interaction(cort,temp))) +
  stat_smooth(method = "glm", span=3, formula = y ~ x, se = TRUE, linewidth=1, alpha = 0.1) +
  scale_color_manual(values = c("CORT.Cold"="darkblue", "Control.Cold"="cyan", "CORT.Hot"="black", "Control.Hot"="grey"),
    labels=c("CORT-Cold (n=11)", "Control-Cold (n=12)", "CORT-Hot (n=11)", "Control-Hot (n=10)")
  ) +
  facet_wrap(deli_associative$group, scales = "free_y", ncol = 1) +
  theme(strip.placement = "outside") +  
  theme(strip.background = element_blank()) +
  labs(y = "Probability of correct choice", x = "Trial") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),
    axis.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 12, family = "sans")
  )
print(plot_deli)
ggsave("./output/figures/fig_deli_asso.png", plot=plot_deli, width = 15, height = 12, units = "cm", dpi = 3000)
 
# Plot guichenoti associative
plot_guich <- ggplot(guich_associative, aes(x = Associative_Trial, y = FC_associative, color=interaction(cort,temp))) +
  stat_smooth(method = "glm", span=3, formula = y ~ x, se = TRUE, linewidth=1, alpha = 0.1) +
  scale_color_manual(values = c("CORT.Cold"="darkblue", "Control.Cold"="cyan", "CORT.Hot"="black", "Control.Hot"="grey"),
    labels=c("CORT-Cold (n=11)", "Control-Cold (n=12)", "CORT-Hot (n=11)", "Control-Hot (n=10)")
  ) +
  facet_wrap(guich_associative$group, scales = "free_y", ncol = 1) +
  theme(strip.placement = "outside") +
  theme(strip.background = element_blank()) +
  labs(y = "Probability of correct choice", x = "Trial") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),
    axis.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 12, family = "sans")
  )
print(plot_guich)
ggsave("./output/figures/fig_guich_asso.png", plot=plot_guich, width = 15, height = 12, units = "cm", dpi = 3000)

# Plot delicata reversal
label_mapping <- c("R_B" = "Blue", "B_R" = "Red")
deli_reversal$group <- factor(deli_reversal$group, levels = levels(deli_reversal$group), labels = label_mapping)
deli_rev <- ggplot(deli_reversal, aes(x = trial_reversal, y = FC_reversal, color=interaction(cort,temp))) +
  stat_smooth(method = "glm", formula = y ~ x, se = TRUE, linewidth=1, alpha = 0.1) +
  scale_color_manual(values = c("CORT.Cold"="darkblue", "Control.Cold"="cyan", "CORT.Hot"="black", "Control.Hot"="grey"),
    labels=c("CORT-Cold (n=11)", "Control-Cold (n=12)", "CORT-Hot (n=11)", "Control-Hot (n=10)")
  ) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  theme(strip.placement = "outside") +  
  theme(strip.background = element_blank()) +
  labs(y="", x = "Trial", color="Treatments") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),
    axis.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 12, family = "sans")
  )
print(deli_rev)
ggsave("./output/figures/fig_deli_rev.png", plot=deli_rev, width = 15, height = 12, units = "cm", dpi = 3000)

# Plot guichenoti reversal
label_mapping <- c("R_B" = "Blue", "B_R" = "Red")
guich_reversal$group <- factor(guich_reversal$group, levels = levels(deli_reversal$group), labels = label_mapping)
guich_rev <- ggplot(guich_reversal, aes(x = trial_reversal, y = FC_reversal, color=interaction(cort,temp))) +
  stat_smooth(method = "glm", formula = y ~ x, se = TRUE, linewidth=1, alpha = 0.1) +
  scale_color_manual(values = c("CORT.Cold"="darkblue", "Control.Cold"="cyan", "CORT.Hot"="black", "Control.Hot"="grey"),
    labels=c("CORT-Cold (n=11)", "Control-Cold (n=12)", "CORT-Hot (n=11)", "Control-Hot (n=10)")
  ) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  theme(strip.placement = "outside") +
  theme(strip.background = element_blank()) +
  labs(y="", x = "Trial", color="Treatments") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),
    axis.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(size = 14, family = "sans"),
    legend.text = element_text(size = 12, family = "sans")
  )
print(guich_rev)
ggsave("./output/figures/fig_guich_rev.png", plot=guich_rev, width = 15, height = 12, units = "cm", dpi = 3000)