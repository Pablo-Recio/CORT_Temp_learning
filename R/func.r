#################### 
# PACKAGES
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes)
#
#################### 
# MODELS
####################
#
# Fit Associative models and extract posteriors both per each treatment:
#' @title fit_asso Function
#' @description Fit brm models for different the associative task
#' @param sp To select the species of interest
#' @param bias To select the 'group' we want, i.e. to analyse those who learn Red first and Blue second
#' @return Posteriors of fitted brm model for each treatment
fit_asso <- function(sp, bias) {
  # Subsample the data
  data <- data_associative %>%
    group_by(lizard_id) %>%
      filter(species %in% sp, group %in% bias) %>%
    ungroup()
  data.frame()
  #Define levels of treatment
  treatment_levels <- c("Control Cold", "CORT Cold", "Control Hot", "CORT Hot")
   # Initialize an empty data frame for result_df
  result_df <- data.frame()
  # Use purrr::map_dfr for iteration and result_df construction
  result_df <- for(treatment_level in treatment_levels) {
    # Subset data for the current treatment level
    data$trt <- relevel(data$trt, ref = treatment_level)
    # Fit the model
    model <- brm(FC_associative ~ Associative_Trial * trt + (1 + Associative_Trial | lizard_id),
                data = data,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
    # Extract Associative_Trial posteriors
    posterior <- as_draws(model)
    posterior_df <- as.data.frame(posterior)
    associative_trial_samples <- posterior_df[, grepl("b_Associative_Trial|b_Intercept", colnames(posterior_df))]
    # Make data frame with all posteriors including treatment_level
    df <- data.frame(associative_trial_samples,
          treatment_level = treatment_level)
    return(df)
    }
return(result_df)
}
#
# Fit Reversal models function:
#' @title fit_rev Function
#' @description Fit brm models for different the reversal task
#' @param sp To select the species of interest
#' @param bias To select the 'group' we want, i.e. to analyse those who learn Red first and Blue second
#' @return Posteriors of fitted brm model for each treatment
fit_rev <- function(sp, bias) {
  # Subsample the data
  data <- data_reversal %>%
    group_by(lizard_id) %>%
      filter(species %in% sp, group %in% bias) %>%
    ungroup()
  data.frame()
  #Define levels of treatment
  treatment_levels <- unique(data$trt)
  # Initialize an empty data frame for result_df
  result_df <- data.frame()
  # Use purrr::map_dfr for iteration and result_df construction
  result_df <- purrr::map_dfr(treatment_levels, function(treatment_level) {
    # Subset data for the current treatment level
    data$trt <- relevel(data$trt, treatment_level)
    # Fit the model
    model <- brm(FC_reversal ~ trial_reversal*trt + (1 + trial_reversal | lizard_id),
                data = data,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
    # Extract Associative_Trial posteriors
    posterior <- as_draws(model)
    posterior_df <- as.data.frame(posterior)
    reversal_trial_samples <- posterior_df[, grepl("b_Intercept", "b_trial_reversal", colnames(posterior_df))]
    # Make data frame with all posteriors including treatment_level
    df <- data.frame(reversal_trial_samples,
          treatment_level = treatment_level)
    return(df)
    })
return(result_df)
}
#
# Extract posteriors function:
#' @title extract_posteriors function
#' @description Obtain posteriors for all factors from different models
#' @param model Path to the model
#' @return Data frame with posteriors
extract_posteriors <- function(model) {
  # Extract posteriors
  posterior <- posterior_samples(model, pars = "^b")
  # Calculate means
  means <- plyr::ldply(lapply(posterior, mean))
  # Calculate quantiles
  quantiles <- plyr::ldply(lapply(posterior, function(x) quantile(x, c(0.025, 0.975))))
  # Combine means, quantiles, and p-values into a single data frame
  result <- cbind(means, quantiles)
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


# Extract the mean of Associative_Trial estimates
# row_means <- rowMeans(associative_trial_samples)
# global_mean <- mean(row_means)
# Extract the quantiles of Associative_Trial estimates
# row_quantiles <- apply(associative_trial_samples, 1, function(x) quantile(x, c(0.025, 0.975)))
# global_quantiles <- quantile(row_quantiles, c(0.025, 0.975))
# Make df data frame
# df <- data.frame(
  #treatment_level = treatment_level,
  #global_mean = global_mean,
  #quantile_0.025 = global_quantiles[1],
  #quantile_0.975 = global_quantiles[2])