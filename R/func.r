#################### 
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x)
#
####################
####################
# Extract sample size per group
#' @title sample
#' @description Estract sample size per group
#' @param sp To select the species of interest ("deli"/"guich")
#' @param bias To select group depending on the colour assigned as correct for each task ("blue"/"red")
#' @param corti To select cort treatment ("CORT"/"Control")
#' @param therm To select temp treatment ("Cold"/"Hot")
sample <- function(sp, bias, corti, therm){
  #Specify database
  data <- data_clean
  # Count sample
  sample_size <- data %>%
                filter(species == sp, group == bias, cort == corti, temp == therm, Trial == 1) %>%
                group_by(lizard_id) %>%
                summarise(n = n()) %>%
                summarise(total_count = sum(n)) %>%
                pull(total_count)
  return(sample_size)
}
####################
####################
# Fit models and extract posteriors both per each treatment:
#' @title fit_m Function
#' @description Fit brm models for the associative task
#' @param sp To select the species of interest ("deli"/"guich")
#' @param bias To select group depending on the colour assigned as correct for each task ("blue"/"red")
#' @param refir To choose whether to refit the models (TRUE, default) or use the ones already made (FALSE)
#' @return Raw posteriors of fitted brm model for each treatment, species, and group (df)
fit_m <- function(sp, bias, refit = TRUE) {
  data <- data_clean
  formula <- FC_associative ~ Trial*cort*temp + (1 + Trial|lizard_id)
  #Specify species
    if (sp == "deli"){
      sp_data <- data %>%
            group_by(lizard_id) %>%
            filter(species == "delicata") %>%
            ungroup() %>%
      data.frame() 
    } else {
      if(sp == "guich"){
        sp_data <- data %>%
              group_by(lizard_id) %>%
              filter(species == "guichenoti") %>%
              ungroup() %>%
        data.frame()
      } else {
        stop("Species not valid")
      }
    }
  #Specify bias/group
    if (bias == "blue"){
      sub_data <- sp_data %>%
            group_by(lizard_id) %>%
            filter(group == "Blue") %>%
            ungroup() %>%
      data.frame() 
    } else {
      if(bias == "red"){
        sub_data <- sp_data %>%
              group_by(lizard_id) %>%
              filter(group == "Red") %>%
              ungroup() %>%
        data.frame()
      } else {
        stop("Group/colour not valid")
      }
    }
  #Fit the model only if it has not been fit yet (if refit=TRUE)
  if(refit){
    # Fit the model
    model <- brm(formula,
                data = sub_data,
                family = bernoulli(link = "logit"),
                chains = 4, cores = 4, iter = 3000, warmup = 1000, control = list(adapt_delta = 0.99))
    # Write the model to a file
    saveRDS(model, file = paste0(here("output/models/"), sp, "_", bias, ".rds"))
  } else {
      # Read the model from a file
      model <- readRDS(file = paste0(here("output/models/"), sp, "_", bias, ".rds"))
  } 
  # Extract posteriors
  posteriors <- as_draws_df(model)
  return(posteriors)
}
###################
###################
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
####################
####################
# Estimate p-values using pmcm
#' @title pMCMC Function
#' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
#' @param null A numeric value decsribing what the null hypothesis should be
#' @param twotail Whether to conduct a one-tailed hypothesis or a two-tailed hypotheses. Default = true indicating a two-tailed test will be done.
#' @param dir The direction of the one-tail test (<) or (>)
pmcmc <- function(x, null = 0, twotail = TRUE, dir){
  if(twotail){
    2*(1 - max(table(x<=null) / length(x)))
  } else{
    if(dir == "<"){
    (max(sum(x>=null) / length(x)))
    } else{
      if(dir == ">"){
        (max(sum(x<=null) / length(x)))
      } else{
        stop("dir not valid")
      }
    }
  }
}
####################
####################
# Function to format numbers with 2 decimal places
#' @title format_dec
#' @param x The object
#' @param n The number of decimals
format_dec <- function(x, n) {
  z <- sprintf(paste0("%.",n,"f"), x)
  return(as.numeric(z))
}
####################
####################
# Function to create each df for the figure
#' @title df_fig
#' @param df to select the df
#' @param specie to add data of the specie
#' @param group to add data of group
df_fig <- function(data_fig, specie, group){
  # Select the original df
  df <- data_fig
  # Create a vector with treatments
  treatments <- c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
  # Create new matrix and new dfs
  treat_matrix <- matrix(NA, nrow = 8000, ncol = 35)
  fig_df <- data.frame()
  # Loop through treatments
  for(t in treatments){
    if(t == "CORT-Cold"){
      m <- df$b_Trial
      u <- df$b_Intercept
    } else if(t == "Control-Cold"){
        m <- (df$'b_Trial:cortControl' + df$b_Trial)
        u <- (df$b_cortControl + df$b_Intercept)
    } else if(t == "CORT-Hot"){
        m <- (df$'b_Trial:tempHot' + df$b_Trial)
        u <- (df$b_tempHot + df$b_Intercept)
    } else if(t == "Control-Hot"){
        m <- (df$'b_Trial:cortControl:tempHot' + df$b_Trial+ df$'b_Trial:cortControl' + df$'b_Trial:tempHot')
        u <- (df$'b_cortControl:tempHot' + df$b_cortControl + df$b_tempHot + df$b_Intercept)
    } else {
    stop("loop wrong")
    }
    # Loop per treatment
    for(x in 0:35){
      for(j in 1:8000){
        value <- exp(u[j] + m[j] * x) / (1 + exp(u[j] + m[j] * x))
        treat_matrix[j, x] <- value
      }
    }
    treat_df <- as.data.frame(treat_matrix)
    colnames(treat_df) <- paste0("X", 1:35)  # Adjust column names
    treat_df <- gather(treat_df, key = "Trial", value = "Value")  # Reshape data frame
    treat_df$Treatment <- t
    fig_df <- rbind(fig_df, treat_df)
  }
  # Add species and group information
  fig_df$Species <- specie
  fig_df$Group <- group
  return(fig_df)
}
####################
####################
# Function to create the plot
#' @title plotting
#' @param df to select the df
plotting <- function(df){
  plot <- ggplot(df, aes(x = Trial, y = Mean_Predicted_prob, color = Treatment)) +
  geom_line(linewidth = 1.75) +
  scale_color_manual(values = c("CORT-Cold"="darkblue", "Control-Cold"="cyan", "CORT-Hot"="black", "Control-Hot"="#616161")) +
  geom_ribbon(aes(ymin = Mean_Predicted_prob - SE_Predicted_prob, ymax = Mean_Predicted_prob + SE_Predicted_prob, fill = Treatment), color = NA, alpha = 0.075) + 
  scale_fill_manual(values = c("CORT-Cold"="darkblue", "Control-Cold"="cyan", "CORT-Hot"="black", "Control-Hot"="#616161")) +
  theme_classic() +
  facet_grid2(Species ~ Group, scale = "free_x", space = "free_x", axes = "all") +
  theme(strip.placement = "outside") +  
  theme(strip.background = element_blank()) +
  labs(y = "Predicted probability of correct choice", x = "Trial") +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm")) + 
  theme(
    axis.title = element_text(size = 15, family = "Times New Roman"),
    axis.text = element_text(size = 10, family = "Times New Roman"),
    legend.title = element_text(size = 12, family = "Times New Roman"),
    legend.text = element_text(size = 10, family = "Times New Roman"),
    strip.text = element_text(size = 13, family = "Times New Roman")
  )  
  return(plot)
}