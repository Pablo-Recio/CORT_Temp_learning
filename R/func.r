#################### 
# PACKAGES
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes)
#
#################### 
# MODELS & RAW POSTERIORS
####################
#
# Fit models and extract posteriors both per each treatment:
#' @title fit_m Function
#' @description Fit brm models for different the associative task
#' @param type To define if we are analysing the reversal or the associative
#' @param sp To select the species of interest
#' @param bias To select the 'group' we want, i.e. to analyse those who learn Red first and Blue second
#' @return Estimates of fitted brm model for each treatment, species, and group (df)
fit_m <- function(type, sp, bias, refit = TRUE) {
  #Specify the type
  if (type == "asso"){
    data <- data_associative
    formula <- FC_associative ~ Associative_Trial + (1 + Associative_Trial|lizard_id)
  }else{
    if(type == "rev"){
      data <- data_reversal
      formula <- FC_reversal ~ trial_reversal + (1 + trial_reversal|lizard_id)
    } else {
      cat("Option not valid\n")
      return(NULL)  # Return NULL in case of an invalid option
    }
  }
  # Create global database for all the posteriors of all the models
  global <- data.frame()
  # Start the loop for all species and groups
  for(species_level in sp) {
    for(bias_level in bias) {
      # Subsample the data
      sub_data <- data %>%
        group_by(lizard_id) %>%
          filter(species_level %in% sp, bias_level %in% bias) %>%
        ungroup()
      data.frame()
      #Define levels of treatment
      treatment_levels <- c("CORT Cold", "Control Cold", "CORT Hot", "Control Hot")
      # Initialize an empty data frame for result_df
      result_df <- data.frame()
      # Use purrr::map_dfr for iteration and result_df construction
      result_df <- purrr::map_dfr(treatment_levels, function(treatment_level) {
        # Subset data for the current treatment level
        sub_data$trt <- treatment_level
        
      if(refit){
        # Fit the model
        model <- brm(formula,
                    data = sub_data,
                    family = bernoulli(link = "logit"),
                    chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
        
        # Write the model to a file'
        saveRDS(model, file = paste0(here("output/models/"), type, "_", species_level, "_", bias_level, "_", treatment_level, ".rds"))
        } else {
          # Read the model from a file
          model <- readRDS(file = paste0(here("output/models/"), type, "_", species_level, "_", bias_level, "_", treatment_level, ".rds"))
        } 
        # Extract estimates of "trial" and the Intercept
        posterior <- as_draws(model)
        posterior_df <- as.data.frame(posterior)
        trial_samples <- posterior_df[, grepl("b_Associative_Trial|b_Intercept", colnames(posterior_df))]
        # Make data frame with all posteriors including treatment_level
        df <- data.frame(trial_samples,
              treatment_level = treatment_level,
              species_level = species_level,
              bias_level = bias_level)
        return(df)
       })
    global <- bind_rows(global, result_df)
    }
  }
return(global)
}
#
# Tidy estimates from fit_m
#' @title tidy_post
#' @description The df obtained by fit_m contains the estimates of the Intercept and fixed effects of "Choice ~ trial (1 + lizard)"
#' for each treatment, species, and group,but it extracts the estimates of the four chains, each chain in a column, here we want to
#' tidy up to two columns (Intercept and trial), plus the ones indicating species, group, tand treatment
#' @param df Dataframe used
#' @return Same df but tidy
tidy_post <- function(df) {
  # Select data
  data <- df
  # Split original df into four by chain, and give the columns a common name
  res1 <- data%>%select(matches("^X1.b_"),treatment_level, species_level,bias_level)%>%rename_all(~sub("^X1.b_", "", .))
  res2 <- data%>%select(matches("^X2.b_"),treatment_level, species_level,bias_level)%>%rename_all(~sub("^X2.b_", "", .))
  res3 <- data%>%select(matches("^X3.b_"),treatment_level, species_level,bias_level)%>%rename_all(~sub("^X3.b_", "", .))
  res4 <- data%>%select(matches("^X4.b_"),treatment_level, species_level,bias_level)%>%rename_all(~sub("^X4.b_", "", .))
  # Bind data again
  new_df <- bind_rows(res1, res2,res3,res4)
return(new_df) 
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