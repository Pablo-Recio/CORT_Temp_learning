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
#| label: setup
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x, PupillometryR, cowplot, png, grid)
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
#| label: fig-Methods
#| fig.cap: "Experimental design of early environment manipulation (**A**) and learning tasks (**B**). Stages 1-3 indicate the different phases of the habituation process. In the associative and reversal tasks, white lids show the ramps where the food reward was not attainable."

knitr::include_graphics("./Others/LEARN_FIG_1.png")

#
#
#
#
#
#
#
#
#| label: cleaning
source(here("R", "1_data_process.R"))
# Once we run this code, the result will be the main data base used (data_clean).
write.csv(data_clean, file= "./output/databases_clean/data_clean.csv")
#
#
#
#| label: sampleSize
# List with the sample sizes from the database (here("output/databases_clean/data_clean.csv") as the sample size per species and group is the same on each task. We used function sample (see func.R) to estimate the sample size per treatment and species.
source(here("R", "func.R"))
#
specie <- c("delicata", "guichenoti")
groups <- c("Red", "Blue")
hormone <- c("CORT", "Control")
temperature <- c("Cold", "Hot")
#
n_list <- list()
#
for(i in 1:length(specie)){
  for(j in 1:length(groups)){
    for(k in 1:length(hormone)){
      for(l in 1:length(temperature)){
        n_list[[paste(specie[i], groups[j], hormone[k], temperature[l], sep = "_")]] <- sample(specie[i], groups[j], hormone[k], temperature[l])
      }
    }
  }
}
#
#
#
#
#| label: models
# Fitting the model and extraction of posteriors for both species using fit_m function (see func.r in R folder). The result everytime the function is used is a df with the posteriors of the model. The functions saves the model automatically in output/models; and when the parameter refit = FALSE then the posteriors are extracted from the model previously written instead of fitting the model again each time.
source(here("R", "func.R"))
# Model formula: FC_associative ~ Trial*cort*temp + (1 + Trial|lizard_id) + (1| clutch)
## A.1) L. delicata
deli_red <- fit_m("deli", "red", refit = FALSE)
write.csv(deli_red, file= "./output/Checking/deli_red.csv")
deli_blue <- fit_m("deli", "blue", refit = FALSE)
write.csv(deli_blue, file= "./output/Checking/deli_blue.csv")
## A.2) L. guichenoti
guich_red <- fit_m("guich", "red", refit = FALSE)
write.csv(guich_red, file= "./output/Checking/guich_red.csv")
guich_blue <- fit_m("guich", "blue", refit = FALSE)
write.csv(guich_blue, file= "./output/Checking/guich_blue.csv")
#
#
#
#| label: results
# Rename some of the posteriors and obtain estimates for the "learning rate" for each group and treatment. This "learning rate" consists on the slope of the posterior of the Trial parameter. The slope is the rate at which the probability of choosing the correct feeder increases with each trial. The slope is obtained by adding the posterior of the Trial parameter to the posterior of the interaction between Trial and the treatment (cort and temp). 
## 1) L. delicata
### Group = red
dar_CORTCold <- deli_red$b_Trial #Slope for treatment CORT-Cold for L. delicata, red as the right choice
dar_ControlCold <- (deli_red$'b_Trial:cortControl' + deli_red$b_Trial)
dar_CORTHot <- (deli_red$'b_Trial:tempHot' + deli_red$b_Trial)
dar_ControlHot <- (deli_red$'b_Trial:cortControl:tempHot' + deli_red$b_Trial+ deli_red$'b_Trial:cortControl' + deli_red$'b_Trial:tempHot')
### Group = blue
dab_CORTCold <- deli_blue$b_Trial
dab_ControlCold <- (deli_blue$'b_Trial:cortControl' + deli_blue$b_Trial)
dab_CORTHot <- (deli_blue$'b_Trial:tempHot' + deli_blue$b_Trial)
dab_ControlHot <- (deli_blue$'b_Trial:cortControl:tempHot' + deli_blue$b_Trial + deli_blue$'b_Trial:cortControl' + deli_blue$'b_Trial:tempHot')
## 2) L. guichenoti
### Group = red
gar_CORTCold <- guich_red$b_Trial
gar_ControlCold <- (guich_red$'b_Trial:cortControl' + guich_red$b_Trial)
gar_CORTHot <- (guich_red$'b_Trial:tempHot' + guich_red$b_Trial)
gar_ControlHot <- (guich_red$'b_Trial:cortControl:tempHot' + guich_red$b_Trial + guich_red$'b_Trial:cortControl' + guich_red$'b_Trial:tempHot')
### Group = blue
gab_CORTCold <- guich_blue$b_Trial
gab_ControlCold <- (guich_blue$'b_Trial:cortControl' + guich_blue$b_Trial)
gab_CORTHot <- (guich_blue$'b_Trial:tempHot' + guich_blue$b_Trial)
gab_ControlHot <- (guich_blue$'b_Trial:cortControl:tempHot' + guich_blue$b_Trial + guich_blue$'b_Trial:cortControl' + guich_blue$'b_Trial:tempHot')
#
#
#
#
#
#
#
#| label: fig-deli
#| fig.cap: "Results for Lampropholis delicata for both colour groups "red" (A,B) and "blue" (C, D). Panels A and C show the predicted probability of choosing the correct feeder first over trials. The lines represent the mean predicted probability of choosing the correct feeder first on each trial, and the shaded areas indicate the standard deviation of the mean; both were obtained by using the slope and intercept estimates from the posterior distributions. The different colours indicate the different treatments. Panels B and D show the distribution of the estimates of slopes per each treatment. The x-axis represents the slope estimate, and in the y-axis are the density of the estimates. The different colours indicate the different treatments. Points and bars represent the mean and standard deviation of the mean of the estimates, respectively."
source(here("R", "func.R"))
# First step, create the dfs for all models
## Red
fig_dar_AC_df <- df_plotAC("deli_red")
fig_dar_BD1_df <- df_plotBD1("deli_red")
fig_dar_BD2_df <- df_plotBD2("deli_red")
## Blue
fig_dab_AC_df <- df_plotAC("deli_blue")
fig_dab_BD1_df <- df_plotBD1("deli_blue")
fig_dab_BD2_df <- df_plotBD2("deli_blue")
#
# Second step, create the plots
fig_dar <- plotting(sp = "deli", col = "red", df_prob = fig_dar_AC_df, df_violin = fig_dar_BD1_df, df_points = fig_dar_BD2_df)
fig_dab <- plotting(sp = "deli", col = "blue", df_prob = fig_dab_AC_df, df_violin = fig_dab_BD1_df, df_points = fig_dab_BD2_df)
#
# Third step, merge plots
fig_deli <- plot_grid(fig_dar, fig_dab, ncol = 1) +
  theme(plot.margin = margin(1, 1, 1.5, 1, "cm")) +
  annotate("text", x = 0.5, y = 0.985, label = "Lampropholis delicata", size = 4.5, color = "black", fontface = "italic")
ggsave("./output/figures/fig_deli.png", plot=fig_deli, width = 25, height = 18, units = "cm", dpi = 600)
knitr::include_graphics("./output/figures/fig_deli.png")
#
#
#
#
#
#| label: fig-guich
#| fig.cap: "Results for Lampropholis guichenoti for both colour groups "red" (A,B) and "blue" (C, D). Panels A and C show the predicted probability of choosing the correct feeder first over trials. The lines represent the mean predicted probability of choosing the correct feeder first on each trial, and the shaded areas indicate the standard deviation of the mean; both were obtained by using the slope and intercept estimates from the posterior distributions. The different colours indicate the different treatments. Panels B and D show the distribution of the estimates of slopes per each treatment. The x-axis represents the slope estimate, and in the y-axis are the density of the estimates. The different colours indicate the different treatments. Points and bars represent the mean and standard deviation of the mean of the estimates, respectively."
source(here("R", "func.R"))
# First step, create the dfs for all models
## Red
fig_gar_AC_df <- df_plotAC("guich_red")
fig_gar_BD1_df <- df_plotBD1("guich_red")
fig_gar_BD2_df <- df_plotBD2("guich_red")
## Blue
fig_gab_AC_df <- df_plotAC("guich_blue")
fig_gab_BD1_df <- df_plotBD1("guich_blue")
fig_gab_BD2_df <- df_plotBD2("guich_blue")
#
# Second step, create the plots
fig_gar <- plotting(sp = "guich", col = "red", df_prob = fig_gar_AC_df, df_violin = fig_gar_BD1_df, df_points = fig_gar_BD2_df)
fig_gab <- plotting(sp = "guich", col = "blue", df_prob = fig_gab_AC_df, df_violin = fig_gab_BD1_df, df_points = fig_gab_BD2_df)
#
# Third step, merge plots
fig_guich <- plot_grid(fig_gar, fig_gab, ncol = 1) +
  theme(plot.margin = margin(1, 1, 1.5, 1, "cm")) +
  annotate("text", x = 0.5, y = 0.985, label = "Lampropholis guichenoti", size = 4.5, color = "black", fontface = "italic")
ggsave("./output/figures/fig_guich.png", plot=fig_guich, width = 25, height = 18, units = "cm", dpi = 600)
knitr::include_graphics("./output/figures/fig_guich.png")
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
cat("\\newpage")
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
#| label: tbl-bias
#| tbl-cap: "Probability of choosing the correct feeder in the first trial when the correct feeder was blue (Prob Blue) or red (Prob Red) for each species and each treatment. p~mcmc~ tested the hypothesis that the probability is >0.33. In bold, those values that are significant (p~mcmc~ <0.05)"
source(here("R", "func.R"))
# First we estimate the probability of choosing right in the first trial using the intercepts from the posteriors
## 1) L. delicata
### Group = red
probright_drCORTCold <- exp(deli_red$b_Intercept)/(1+exp(deli_red$b_Intercept))
probright_drControlCold <- exp(deli_red$b_cortControl + deli_red$b_Intercept)/(1+exp(deli_red$b_cortControl + deli_red$b_Intercept))
probright_drCORTHot <- exp(deli_red$b_tempHot + deli_red$b_Intercept)/(1+exp(deli_red$b_tempHot + deli_red$b_Intercept))
probright_drControlHot <- exp(deli_red$'b_cortControl:tempHot' + deli_red$b_cortControl + deli_red$b_tempHot + deli_red$b_Intercept)/(1+exp(deli_red$'b_cortControl:tempHot' + deli_red$b_cortControl + deli_red$b_tempHot + deli_red$b_Intercept))
### Group = blue
probright_dbCORTCold <- exp(deli_blue$b_Intercept)/(1+exp(deli_blue$b_Intercept))
probright_dbControlCold <- exp(deli_blue$b_cortControl + deli_blue$b_Intercept)/(1+exp(deli_blue$b_cortControl + deli_blue$b_Intercept))
probright_dbCORTHot <- exp(deli_blue$b_tempHot + deli_blue$b_Intercept)/(1+exp(deli_blue$b_tempHot + deli_blue$b_Intercept))
probright_dbControlHot <- exp(deli_blue$'b_cortControl:tempHot' + deli_blue$b_cortControl + deli_blue$b_tempHot + deli_blue$b_Intercept)/(1+exp(deli_blue$'b_cortControl:tempHot' + deli_blue$b_cortControl + deli_blue$b_tempHot + deli_blue$b_Intercept))
## 2) L. guichenoti
### Group = red
probright_grCORTCold <- exp(guich_red$b_Intercept)/(1+exp(guich_red$b_Intercept))
probright_grControlCold <- exp(guich_red$b_cortControl + guich_red$b_Intercept)/(1+exp(guich_red$b_cortControl + guich_red$b_Intercept))
probright_grCORTHot <- exp(guich_red$b_tempHot + guich_red$b_Intercept)/(1+exp(guich_red$b_tempHot + guich_red$b_Intercept))
probright_grControlHot <- exp(guich_red$'b_cortControl:tempHot' + guich_red$b_cortControl + guich_red$b_tempHot + guich_red$b_Intercept)/(1+exp(guich_red$'b_cortControl:tempHot' + guich_red$b_cortControl + guich_red$b_tempHot + guich_red$b_Intercept))
### Group = blue
probright_gbCORTCold <- exp(guich_blue$b_Intercept)/(1+exp(guich_blue$b_Intercept))
probright_gbControlCold <- exp(guich_blue$b_cortControl + guich_blue$b_Intercept)/(1+exp(guich_blue$b_cortControl + guich_blue$b_Intercept))
probright_gbCORTHot <- exp(guich_blue$b_tempHot + guich_blue$b_Intercept)/(1+exp(guich_blue$b_tempHot + guich_blue$b_Intercept))
probright_gbControlHot <- exp(guich_blue$'b_cortControl:tempHot' + guich_blue$b_cortControl + guich_blue$b_tempHot + guich_blue$b_Intercept)/(1+exp(guich_blue$'b_cortControl:tempHot' + guich_blue$b_cortControl + guich_blue$b_tempHot + guich_blue$b_Intercept))
#
# Second we build a df with all the mean probabilities of choosing right in the first trial for each treatment and species
prob_df <- data.frame(
  probright_drCORTCold_t = format_dec(mean(probright_drCORTCold), 3),
  probright_drControlCold_t = format_dec(mean(probright_drControlCold), 3),
  probright_drCORTHot_t = format_dec(mean(probright_drCORTHot), 3),
  probright_drControlHot_t = format_dec(mean(probright_drControlHot), 3),
  probright_dbCORTCold_t = format_dec(mean(probright_dbCORTCold), 3),
  probright_dbControlCold_t = format_dec(mean(probright_dbControlCold), 3),
  probright_dbCORTHot_t = format_dec(mean(probright_dbCORTHot), 3),
  probright_dbControlHot_t = format_dec(mean(probright_dbControlHot), 3),
  probright_grCORTCold_t = format_dec(mean(probright_grCORTCold), 3),
  probright_grControlCold_t = format_dec(mean(probright_grControlCold), 3),
  probright_grCORTHot_t = format_dec(mean(probright_grCORTHot), 3),
  probright_grControlHot_t = format_dec(mean(probright_grControlHot), 3),
  probright_gbCORTCold_t = format_dec(mean(probright_gbCORTCold), 3),
  probright_gbControlCold_t = format_dec(mean(probright_gbControlCold), 3),
  probright_gbCORTHot_t = format_dec(mean(probright_gbCORTHot), 3),
  probright_gbControlHot_t = format_dec(mean(probright_gbControlHot), 3)) %>%
gather(key = "column", value = "value", .) %>%
rename("prob_choice" = "value") %>%
data.frame()
write.csv(prob_df, file= "./output/Checking/prob_df.csv")
# Then we build a df with the statistical test (one-tailed pmcmc) that the probability is >0.33 (the probability of choosing right by chance)
comp_pval <- data.frame(
  probright_drCORTCold_t = format_p(pmcmc(probright_drCORTCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_drControlCold_t = format_p(pmcmc(probright_drControlCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_drCORTHot_t = format_p(pmcmc(probright_drCORTHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_drControlHot_t = format_p(pmcmc(probright_drControlHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_dbCORTCold_t = format_p(pmcmc(probright_dbCORTCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_dbControlCold_t = format_p(pmcmc(probright_dbControlCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_dbCORTHot_t = format_p(pmcmc(probright_dbCORTHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_dbControlHot_t = format_p(pmcmc(probright_dbControlHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_grCORTCold_t = format_p(pmcmc(probright_grCORTCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_grControlCold_t = format_p(pmcmc(probright_grControlCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_grCORTHot_t = format_p(pmcmc(probright_grCORTHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_grControlHot_t = format_p(pmcmc(probright_grControlHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_gbCORTCold_t = format_p(pmcmc(probright_gbCORTCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_gbControlCold_t = format_p(pmcmc(probright_gbControlCold, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_gbCORTHot_t = format_p(pmcmc(probright_gbCORTHot, null = 0.33, twotail = FALSE, dir=">"), 3),
  probright_gbControlHot_t = format_p(pmcmc(probright_gbControlHot, null = 0.33, twotail = FALSE, dir=">"), 3)) %>%
gather(key = "column", value = "value", .) %>%
rename("p_value" = "value") %>%
data.frame()
write.csv(comp_pval, file= "./output/Checking/comp_pval.csv")
# Now we merge them and organise them to have the final dt for the table
table_df <- merge(prob_df, comp_pval, by ="column") %>%
  mutate(Specie = gsub(".*_(d|g).*", "\\1",column)) %>%
  mutate(Group = gsub(".*_(dr|db|gr|gb).*", "\\1",column)) %>%
  mutate(Treatment = gsub(".*(CORTCold|ControlCold|CORTHot|ControlHot).*", "\\1",column)) %>% #Last three commands used to get the species, group and treatment from the column names
  mutate(Treatment = factor(Treatment, 
                    levels = c("CORTCold", "ControlCold", "CORTHot", "ControlHot"),
                    labels = c("CORTCold" = "CORT-Cold", "ControlCold" = "Control-Cold", "CORTHot" = "CORT-Hot", "ControlHot" = "Control-Hot"))) %>%
  mutate(Specie = factor(Specie,
                levels = c("d", "g"),
                labels = c("d" = "L. delicata", "g" = "L. guichenoti"))) %>%
  mutate(Group = recode(Group, "dr" = "Red", "db" = "Blue", "gr" = "Red", "gb" = "Blue")) %>%
  select(-column) %>%
  pivot_wider(names_from = Group, values_from = c(prob_choice, p_value)) %>% #To split the df in Red and Blue
  select(Specie, Treatment, prob_choice_Blue, p_value_Blue, prob_choice_Red, p_value_Red) %>%
  arrange(Specie, Treatment) %>% # To arrange the rows the way I wantdata.frame()
data.frame()
write.csv(table_df, file= "./output/Checking/table_df.csv")
#################
colour_table <- flextable(table_df) %>%
  bold(~ `p_value_Blue` < 0.05, ~ `p_value_Blue` + prob_choice_Blue) %>%
  bold(~ `p_value_Red` < 0.05, ~ `p_value_Red` + prob_choice_Red) %>%
  set_table_properties(width = 1) %>%
  align(align="center", part="all") %>% 
  set_header_labels(prob_choice_Blue = "Prob Blue", 
                    p_value_Blue = "p-value Blue", 
                    prob_choice_Red = "Prob Red", 
                    p_value_Red = "p-value Red") %>%
  italic(j = 1, italic = TRUE, part = "body") %>% # To have names of species in italics
  flextable::compose(i = c(2:4,6:8), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  hline(i = c(4), j = c(1:6), part = "body") %>% # To make some horizontal lines
  vline(j = c(2,4), part = "all") %>% # To make some vertical lines on body
autofit()
colour_table
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#| label: fig-spectrum
#| fig.cap: "Light spectrum of the ramps used in the associative task and the white ones used during habituation. The different colours represent the different ramps."
#
# Load required packages
pacman::p_load(pavo,plotrix,plyr,MCMCglmm,Hmisc,DIZutils)
#
# For importing the files, change the extension of the .jaz files to .txt files!  
# Import spec data for pablos feeders - CHNAGE TO WHERE YOUR FILES ARE LOCATED
specs <- getspec("./Others/Feeders", lim = c(300,700), ext = 'txt') %>%
  procspec(specs, opt = "smooth" , span = 0.2, fixneg = "addmin") 
# Plot smoothed raw data as reflectance curve with approximate rgb colours
plot(specs, type = "o", col = spec2rgb(specs))
#
#
#
#| label: fig-perceived1
#| fig.cap: "Perceived chromatic contrasts between feeders."
#
## Processing the data with the lizard visual model 
# Processing using lizard vision model provided for Dan and Martin's Phyrnocephalus mystaceus paper
liz_vis <- sensmodel(c(360, 440, 493, 571)) 
#
# Apply this over the feeder data visual as the visual system, achromatic is the receptor sensitivity, illum is the illuminate (here using an ideal white illuminate preset), bkg is the background against which we are measuring (here using an ideal background with no influence). Other conditions are applied as used in the P mystaceus paper.
# Does produce an error about luminance but this is noted on the original P mystaceus paper code and should be okay
# summary(feeders_lizardview) can be used to check what conditions were applied in the visual model
feeders_lizardview <- vismodel(specs, visual = liz_vis, achromatic = "l", illum = "ideal", bkg = "ideal", vonkries = TRUE,   relative = FALSE, qcatch = "fi")
#                             
# Calculate colour distances between all feeder measurements
coldistance_feeders <- coldist(feeders_lizardview,achro = TRUE, noise = "neural", n = c(1, 1, 3.5, 6),  weber = 0.1)
#
# Name all contrasts e.g white/green to categorise them
annotated_coldistance_feeders <- coldistance_feeders %>%
  separate(patch1, into = c("color1", "number1"), sep = "_") %>%
  separate(patch2, into = c("color2", "number2"), sep = "_") %>%
  mutate(color_comparison = paste0(color1, "/", color2)) %>%
  mutate(color_comparison = paste0(pmin(color1, color2), "/", pmax(color1, color2))) %>% # Remove duplicate categories(e.g white/green is the same as green/white and all are called green/white)
  subset(!(color_comparison %in% c("B/B", "G/G", "R/R", "W/W")))
#
# Plot the colour differences between the feeders: 
# Note - above 1 is expected to be visually different under the selected visual system (here, lizard view)
# dS, or chromatic contrast
dScomparisons <- ggplot(annotated_coldistance_feeders, aes(x = color_comparison, y = dS, fill = color_comparison)) + 
  geom_boxplot() + 
  theme_classic() +
  ylab("Chromatic Contrast (dS)") +
  xlab("Feeder Colours Comparison") +
 scale_fill_manual(values =c("B/G" = 'steelblue',  "B/R" = 'purple', "B/W" = 'lightblue', "G/R" = 'brown', "G/W" = 'lightgreen', "R/W" = 'red')) +
 theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm")) + 
  theme(
    axis.title = element_text(size = 15, family = "Times New Roman"),
    axis.text = element_text(size = 10, family = "Times New Roman"),
    legend.title = element_text(size = 12, family = "Times New Roman"),
    legend.text = element_text(size = 10, family = "Times New Roman"),
    strip.text = element_text(size = 13, family = "Times New Roman")
  ) 
print(dScomparisons)
#
#
#
#| label: fig-perceived2
#| fig.cap: "Perceived achromatic contrasts between feeders."
#dL, or achromatic contrast
dLcomparisons <- ggplot(annotated_coldistance_feeders, aes(x = color_comparison, y = dL, fill = color_comparison)) + 
  geom_boxplot() + 
  theme_classic() +
  ylab("Achromatic Contrast (dL") +
  xlab("Feeder Colours Comparison") +
 scale_fill_manual(values =c("B/G" = 'steelblue',  "B/R" = 'purple', "B/W" = 'lightblue', "G/R" = 'brown', "G/W" = 'lightgreen', "R/W" = 'red')) +
 theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm")) + 
  theme(
    axis.title = element_text(size = 15, family = "Times New Roman"),
    axis.text = element_text(size = 10, family = "Times New Roman"),
    legend.title = element_text(size = 12, family = "Times New Roman"),
    legend.text = element_text(size = 10, family = "Times New Roman"),
    strip.text = element_text(size = 13, family = "Times New Roman")
  )
print(dLcomparisons)
#
#
#
#
cat("\\newpage")
#
#
#
#
# Chunk for calling all the models and making the residuals
## L. delicata
### Red
mod_dar <- readRDS("./output/models/deli_red.rds")
resid_dar <- residuals(mod_dar)
### Blue
mod_dab <- readRDS("./output/models/deli_blue.rds")
resid_dab <- residuals(mod_dab)
## L. guichenoti
### Red
mod_gar <- readRDS("./output/models/guich_red.rds")
resid_gar <- residuals(mod_gar)
### Blue
mod_gab <- readRDS("./output/models/guich_blue.rds")
resid_gab <- residuals(mod_gab)
#
#
#
#
#
#
#
#
bayes_R2(mod_dar)
plot(mod_dar)
#
#
#
cat("\\newpage")
#
#
#
#
bayes_R2(mod_dab)
plot(mod_dab)
#
#
#
cat("\\newpage")
#
#
#
#
#
bayes_R2(mod_gar)
plot(mod_gar)
#
#
#
cat("\\newpage")
#
#
#
#
bayes_R2(mod_gab)
plot(mod_gab)
#
#
#
cat("\\newpage")
#
#
#
#
#
#| label: modelsage
# Fitting the models for L. delicata and L. guichenoti including age as a covariate
#
data_age <- data_clean %>%
  mutate(age.start = as.numeric(age.start)-mean(as.numeric(age.start)))
refit2 = FALSE
formula <- FC_associative ~ age.start + trial_associative*cort*temp + (1 + trial_associative|lizard_id) + (1|clutch)
species <- c("delicata", "guichenoti")
group <- c("Red", "Blue")
#
models <- list() # Create an empty list to store the models
for(s in species){
    # Subset the data for the species
    data_species <- data_age %>%
    filter(species == s)
  for(g in group){
  # Subset the data for the group
  data_group <- data_species %>%
  filter(group == g)
  # Define model file name
  model_file <- paste0("output/models/model_age_", s, "_", g,".rds")
# Fit the model only if it has not been fit yet (if refit2 = TRUE)
  if(refit2 == TRUE){
    # Fit the model
    model <- brm(formula,
                data = data_group,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 3000, warmup = 1000, control = list(adapt_delta = 0.99))
    # Write the model to a file
    saveRDS(model, file = model_file)
  } else {
      # Read the model from a file
      model <- readRDS(file = model_file)
  }
  # Save the model to the list
  models[[paste0(s, "_", g)]] <- model
  }
}
#
model_age_deli_red <- models[["delicata_Red"]]
model_age_deli_blue <- models[["delicata_Blue"]]
model_age_guich_red <- models[["guichenoti_Red"]]
model_age_guich_blue <- models[["guichenoti_Blue"]]
#
#
#
#
#| label: tbl-agedeli_red
#| tbl-cap: "Summary model_age_deli_red"
deli_age_red <- summary(model_age_deli_red)
# Extract fixed effects names and their corresponding estimates, CIs, etc.
fixed_effects <- deli_age_red$fixed %>%
  mutate_all(~ format_dec(., 2))
Predictors <- rownames(fixed_effects)
deli_age_red_tbl_df <- cbind(Predictors, fixed_effects)
# Convert the fixed effects data frame to a flextable
flextable(deli_age_red_tbl_df) %>%
  set_table_properties(., width = 1) %>%
  fontsize(size = 10, part = "all")
#
#
#
#
#| label: tbl-agedeli_blue
#| tbl-cap: "Summary model_age_deli_blue"
deli_age_blue <- summary(model_age_deli_blue)
# Extract fixed effects names and their corresponding estimates, CIs, etc.
fixed_effects <- deli_age_blue$fixed %>%
  mutate_all(~ format_dec(., 2))
Predictors <- rownames(fixed_effects)
deli_age_blue_tbl_df <- cbind(Predictors, fixed_effects)
# Convert the fixed effects data frame to a flextable
flextable(deli_age_blue_tbl_df) %>%
  set_table_properties(., width = 1) %>%
  fontsize(size = 10, part = "all")
#
#
#
#
#| label: tbl-ageguich_red
#| tbl-cap: "Summary model_age_guich_red"
guich_age_red <- summary(model_age_guich_red)
# Extract fixed effects names and their corresponding estimates, CIs, etc.
fixed_effects <- guich_age_red$fixed %>%
  mutate_all(~ format_dec(., 2))
Predictors <- rownames(fixed_effects)
guich_age_red_tbl_df <- cbind(Predictors, fixed_effects)
# Convert the fixed effects data frame to a flextable
flextable(guich_age_red_tbl_df) %>%
  set_table_properties(., width = 1) %>%
  fontsize(size = 10, part = "all")
#
#
#
#
#| label: tbl-ageguich_blue
#| tbl-cap: "Summary model_age_guich_blue"
guich_age_blue <- summary(model_age_guich_blue)
# Extract fixed effects names and their corresponding estimates, CIs, etc.
fixed_effects <- guich_age_blue$fixed %>%
  mutate_all(~ format_dec(., 2))
Predictors <- rownames(fixed_effects)
guich_age_blue_tbl_df <- cbind(Predictors, fixed_effects)
# Convert the fixed effects data frame to a flextable
flextable(guich_age_blue_tbl_df) %>%
  set_table_properties(., width = 1) %>%
  fontsize(size = 10, part = "all")
#
#
#
#
#