

# Packages
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes)


```{r, models1}
#| label: models1
#| echo: false
#| warning: false
library(emmeans)
library(lme4)
# Generalized Mixed Models delicata
GMM_deli <- glmer(FC_associative ~ Associative_Trial*temp*cort + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = deli_associative, family = binomial(link = "logit"))
summary(GMM_deli)
emmeans(GMM_deli, pairwise ~ temp*Associative_Trial)
emmeans(GMM_deli, pairwise ~ cort*temp*Associative_Trial)

# Generalized Mixed Models guichenoti
GMM_guich <- glmer(FC_associative ~ Associative_Trial*temp*cort + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = guich_associative, family = binomial(link = "logit"), control = glmerControl(maxit = 1000))
summary(GMM_guich)

emmeans(GMM_guich, pairwise ~ temp*cort*Associative_Trial)

# Generalized Mixed Models delicata REVERSAL
GMM_deli_rev <- glmer(FC_reversal ~ trial_reversal*temp*cort + group*trial_reversal + (1 + trial_reversal|lizard_id), data = deli_reversal, family = binomial(link = "logit"))
summary(GMM_deli)

emmeans(GMM_deli_rev, pairwise ~ cort*temp*trial_reversal)
emmeans(GMM_deli_rev, pairwise ~ temp*trial_reversal)

# Generalized Mixed Models guichenoti REVERSAL
MM_guich <- glmer(FC_reversal ~ trial_reversal*temp*cort + group*trial_reversal + (1 + trial_reversal|lizard_id), data = guich_reversal, family = binomial(link = "logit"))
summary(GMM_guich)

emmeans(GMM_guich, pairwise ~ temp*cort*Associative_Trial)
```

```{r, models2}
```

```{r, modelsasso}
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

contrasts <- emmeans(deli_mod, specs = ~ temp:Associative_Trial)
contrast(contrasts, method = "pairwise", by = NULL)
```

```{r, modelsrev}
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
posterior <- posterior_samples(deli_mod, pars = "^b")
plyr::ldply(lapply(posterior, mean))
plyr::ldply(lapply(posterior, function(x) quantile(x, c(0.025, 0.975))))

print(posterior)
temp23 <- c(posterior[,"b_Intercept"], (posterior[,"b_Intercept"]+ posterior[,"b_trtB_23"]))
temp28 <- c((posterior[,"b_Intercept"] + posterior[,"b_trtB_28"]), posterior[,"b_trtB_23"])

pmcmc(temp23 - temp28)
dharma(GMM_deli2)


```

```{r, plots}
#| label: plots
#| echo: false
#| warning: false
 
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

```


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

## Merge


#' @title pMCMC Function
 #' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
 #' @param null A numeric value decsribing what the null hypothesis should be
 #' @param twotail Whether to conduct a one-tailed hypothesis or a two-tailed hypotheses. Default = true indicating a two-tailed test will be done.
pmcmc <- function(x, null = 0, twotail = TRUE){
  if(twotail){
    2*(1 - max(table(x<=null) / length(x)))
  } else{
    (1 - max(table(x<=null) / length(x)))
  }
}





We calculated $R^2$ as follows

$$
R^2 = \frac{SS_{reg}}{SS_{tot}}
$$ {#eq-r2}

@eq-r2


We have a total of `r length(unique(data$Lizard_id))`

```{r,tb1}
#| label: tbl-tb1
#| tbl-cap: "Summary of the data"

tab <- data %>% 
  group_by(Lizard_id) %>% 
  summarise(
    n = n()
  ) 
flextable(tab[2:3,])
```