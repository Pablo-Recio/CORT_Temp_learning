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
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes)
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
#| label: fig-fig1
#| fig-cap: "Experimental design of early environment manipulation (**A**) and learning tasks (**B**). Stages 1-3 indicate the different phases of the habituation process. In the associative and reversal tasks, white lids indicate the ramps where the food reward was not attainable."
knitr::include_graphics("./Others/LEARN_FIG_1.svg") 
#
#
#
#
#
#
#
#| label: cleaning
source(here("R", "1_data_process.R"))
# Once we run this code, the result will be two different dataframes, one for the Associative learning task (data_asso) and the other for the Reversal learning task (data_rev).
write.csv(data_asso, file= "./output/databases_clean/data_asso.csv")
write.csv(data_rev, file= "./output/databases_clean/data_rev.csv")
#
#
#
#| label: models
# fitting the model and extraction of posteriors for both types of task and species using fit_m function
source(here("R", "func.R"))
# A) Associative task. Model formula: FC_associative ~ Associative_Trial*cort*temp + (1 + Associative_Trial|lizard_id)
## A.1) L. delicata
deli_asso_red <- fit_m("asso", "deli", "red")
write.csv(deli_asso_red, file= "./output/Checking/deli_asso_red.csv")
deli_asso_blue <- fit_m("asso", "deli", "blue")
write.csv(deli_asso_blue, file= "./output/Checking/deli_asso_blue.csv")
## A.2) L. guichenoti
guich_asso_red <- fit_m("asso", "guich", "red")
write.csv(guich_asso_red, file= "./output/Checking/guich_asso_red.csv")
guich_asso_blue <- fit_m("asso", "guich", "blue")
write.csv(guich_asso_blue, file= "./output/Checking/guich_asso_blue.csv")
# B) Reversal task. Model formula: FC_reversal ~ trial_reversal*cort*temp + (1 + trial_reversal|lizard_id)
## B.1) L. delicata
deli_rev_red <- fit_m("rev", "deli", "red")
write.csv(deli_rev_red, file= "./output/Checking/deli_rev_red.csv")
deli_rev_blue <- fit_m("rev", "deli", "blue")
write.csv(deli_rev_blue, file= "./output/Checking/deli_rev_blue.csv")
## B.2) L. guichenoti
guich_rev_red <- fit_m("rev", "guich", "red")
write.csv(guich_rev_red, file= "./output/Checking/guich_rev_red.csv")
guich_rev_blue <- fit_m("rev", "guich", "blue")
write.csv(guich_rev_blue, file= "./output/Checking/guich_rev_blue.csv")
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
