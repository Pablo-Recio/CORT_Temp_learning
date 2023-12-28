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
#| label: cleaning
source(here("R", "1_data_process.R"))
#
#
#
#| label: models
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
#
#| label: posteriors
#### The function here is to extract the posteriors of the models fitted in the previous chunk
source(here("R", "func.R"))
posteriors_deli_asso <- extract_posteriors(deli_asso)
emm_results <- emmeans(deli_asso, specs = ~trial_associative:temp:cort)
emm_results <- emmeans(deli_asso, specs = ~trial_associative:group)
emm_results <- emmeans(deli_asso, specs = ~trial_associative:temp)
emm_results <- emmeans(deli_asso, specs = ~trial_associative:cort)
summary(emm_results, infer = TRUE)
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
