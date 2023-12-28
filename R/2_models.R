####################################
# 2_MODELS
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom)
#
######## 2.A) BUILD ALL MODELS
#### The function here is to fit all the (brm) models for future analyses using the databases from output/clean_databases already processed using the code in "1_data_process.R"
source(here("R", "func.R"))
# A) Associative task
### L. delicata
deli_asso<-fit_asso(here("output/databases_clean/deli_associative.csv"))
### L. guichenoti
guich_asso<-fit_asso(here("output/databases_clean/guich_associative.csv"))
# B) Reversal task
### L. delicata
deli_rev<-fit_rev(here("output/databases_clean/deli_reversal.csv"))
### L. guichenoti
guich_rev<-fit_rev(here("output/databases_clean/guich_reversal.csv"))
#
#
######## 2.B) EXTRACT POSTERIORS FOR ALL MODELS
source(here("R", "func.R"))
posteriors_deli_asso <- extract_posteriors(deli_asso)
posteriors_guich_asso <- extract_posteriors(guich_asso)
posteriors_deli_rev <- extract_posteriors(deli_rev)
posteriors_guich_rev <- extract_posteriors(guich_rev)
#
emtrends_results <- emtrends(deli_asso, ~ cort * temp * Associative_Trial|group, var = "Associative_Trial", infer = TRUE)
# Display the results
summary(emtrends_results)

posterior <- posterior_samples(deli_mod, pars = "^b")
plyr::ldply(lapply(posterior, mean))
plyr::ldply(lapply(posterior, function(x) quantile(x, c(0.025, 0.975))))

######## 2.C) pmcmc FOR ALL MODELS

# pmcmc for differences between 'Groups'

# Also pmcmc for differences between both species


######## 2.D) MAKE TABLE WITH POSTERIORS AND pmcmc


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