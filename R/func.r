#################### 
# PACKAGES
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, pmcmc)
#
#################### 
# MODELS
####################
#
# Fit Associative models function:
#' @title fit_asso Function
#' @description Fit brm models for different the associative task
#' @param data_asso Path to the database file
#' @return Fitted brm model
fit_asso <- function(data_asso) {
  # Load the data
  data <- read.csv(data_asso)
  # Fit the model
  model <- brm(FC_associative ~ cort*temp*trial_associative + group*trial_associative + (1 + trial_associative|lizard_id),
              data = data,
              family = bernoulli(link = "logit"),
              chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
return(model)
}
#
# Fit Reversal models function:
#' @title fit_rev Function
#' @description Fit brm models for different the reversal task
#' @param data_rev Path to the database file
#' @return Fitted brm model
fit_rev <- function(data_rev) {
  # Load the data
  data <- read.csv(data_rev)
  # Fit the model
  model <- brm(FC_reversal ~ cort*temp*trial_reversal + group*trial_reversal (1 + trial_reversal|lizard_id),
              data = data,
              family = bernoulli(link = "logit"),
              chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
return(model)
}
#
# Extract posteriors functions:
#' @title extract_posteriors function
#' @description Obtain posteriors for all factors from different models
#' @param model Path to the model
#' @return Data frame with posteriors
posteriors <- function(model) {
  # Extract posteriors
  extract_posteriors <- posterior_samples(model, pars = "^b")
  # Calculate means
  means <- plyr::ldply(lapply(posterior, mean))
  # Calculate quantiles
  quantiles <- plyr::ldply(lapply(posterior, function(x) quantile(x, c(0.025, 0.975))))
  # Calculate p-values using pmcmc
  p_values <- pmcmc(posterior, null = 0, twotail = TRUE)
  # Combine means, quantiles, and p-values into a single data frame
  result <- cbind(means, quantiles, p_values)
  return(result)
}


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
