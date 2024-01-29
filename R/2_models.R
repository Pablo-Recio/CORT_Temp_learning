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





```{r, FigResults}
|# label: FigResults
|# fig.cap:
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

```