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
#' @param refit To choose whether to refit the models (TRUE, default) or use the ones already made (FALSE)
#' @return Raw posteriors of fitted brm model for each treatment, species, and group (df)
fit_m <- function(sp, bias, refit = TRUE) {
  data <- data_clean
  formula <- FC_associative ~ Trial*cort*temp + (1 + Trial|lizard_id) + (1|clutch)
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
  return(z)
}
####################
####################
# Function to format p_values
#' @title format_p
#' @param x The object
#' @param n The number of decimals
format_p <- function(x, n) {
  if(equal == TRUE){
    e <- "= "
  } else if(equal == FALSE){
    e <- ""  
  }
  z <- sprintf(paste0("%.",n,"f"), x)
  tmp <- ifelse(as.numeric(z) <= 0.001, "< 0.001",
         ifelse(as.numeric(z) <= 0.05 & as.numeric(z) > 0.001, "< 0.05",
                paste0(e, format_dec(as.numeric(z), 2))))
  return(tmp)
}
####################
####################
# Function to create each df for the probability figure (plots C, F fig-deli and fig-guich)
#' @title df_plotprob
#' @param data to select the df, levels are: "deli_red", "deli_blue", "guich_red", "guich_blue" 
df_plotprob <- function(data){
  if (data == "deli_red"){
    df <- as.data.frame(deli_red)
  } else if (data == "deli_blue"){
    df <- as.data.frame(deli_blue)
  } else if(data == "guich_red"){
    df <- as.data.frame(guich_red)
  } else if(data == "guich_blue"){
    df <- as.data.frame(guich_blue)
  } else {
    stop("Data not valid")
  }
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
  # Transform the database for better use
  fig_df <- fig_df %>%
  mutate(Trial = gsub("X", "", Trial)) %>%
  mutate(Trial = as.numeric(Trial)) %>%
  group_by(Trial, Treatment) %>%
  summarize(
    Mean_Predicted_prob = mean(Value),
    SE_Predicted_prob = sd(Value)
    ) %>%
  ungroup() %>%
data.frame()
  if(data == "deli_red"){
    final_fig_df <- fig_df %>%
      mutate(Treatment = factor(Treatment, 
        levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
        labels = c("CORT-Cold" = paste0("CORT-Cold (n = ", n_list$delicata_Red_CORT_Cold, ")"),
                  "Control-Cold" = paste0("Control-Cold (n = ", n_list$delicata_Red_Control_Cold, ")"),
                  "CORT-Hot" = paste0("CORT-Hot (n = ", n_list$delicata_Red_CORT_Hot, ")"),
                  "Control-Hot" = paste0("Control-Hot (n = ", n_list$delicata_Red_Control_Hot, ")"))
      )) %>%
    data.frame()
  } else if(data == "deli_blue"){
    final_fig_df <- fig_df %>%
      mutate(Treatment = factor(Treatment, 
        levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
        labels = c("CORT-Cold" = paste0("CORT-Cold (n = ", n_list$delicata_Blue_CORT_Cold, ")"),
                  "Control-Cold" = paste0("Control-Cold (n = ", n_list$delicata_Blue_Control_Cold, ")"),
                  "CORT-Hot" = paste0("CORT-Hot (n = ", n_list$delicata_Blue_CORT_Hot, ")"),
                  "Control-Hot" = paste0("Control-Hot (n = ", n_list$delicata_Blue_Control_Hot, ")"))
      )) %>%
    data.frame()
  } else if(data == "guich_red"){
    final_fig_df <- fig_df %>%
      mutate(Treatment = factor(Treatment, 
        levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
        labels = c("CORT-Cold" = paste0("CORT-Cold (n = ", n_list$guichenoti_Red_CORT_Cold, ")"),
                  "Control-Cold" = paste0("Control-Cold (n = ", n_list$guichenoti_Red_Control_Cold, ")"),
                  "CORT-Hot" = paste0("CORT-Hot (n = ", n_list$guichenoti_Red_CORT_Hot, ")"),
                  "Control-Hot" = paste0("Control-Hot (n = ", n_list$guichenoti_Red_Control_Hot, ")"))
      )) %>%
    data.frame()
  } else if(data == "guich_blue"){
    final_fig_df <- fig_df %>%
      mutate(Treatment = factor(Treatment, 
        levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
        labels = c("CORT-Cold" = paste0("CORT-Cold (n = ", n_list$guichenoti_Blue_CORT_Cold, ")"),
                  "Control-Cold" = paste0("Control-Cold (n = ", n_list$guichenoti_Blue_Control_Cold, ")"),
                  "CORT-Hot" = paste0("CORT-Hot (n = ", n_list$guichenoti_Blue_CORT_Hot, ")"),
                  "Control-Hot" = paste0("Control-Hot (n = ", n_list$guichenoti_Blue_Control_Hot, ")"))
      )) %>%
    data.frame()
  }
  return(final_fig_df)
}
####################
####################
# Function to create each df for the slopes and intercept figures (part of plots A, B, D, E fig-deli or fig-guich)
#' @title df_plotviol
#' @param data to select the df, levels are: "deli_red", "deli_blue", "guich_red", "guich_blue"
#' @param var to select whether we use the intercept or the slope (choice/slope)
df_plotviol <- function(data, var){
  if (data == "deli_red"){
      if(var == "choice"){
        lst <- list('CORT-Cold' = probright_drCORTCold,
                  'Control-Cold' = probright_drControlCold,
                  'CORT-Hot' = probright_drCORTHot,
                  'Control-Hot' = probright_drControlHot)
      } else if(var == "slope"){
        lst <- list('CORT-Cold' = dar_CORTCold,
                  'Control-Cold' = dar_ControlCold,
                  'CORT-Hot' = dar_CORTHot,
                  'Control-Hot' = dar_ControlHot)
      }
      df_plotBD <- data.frame(Value = unlist(lst),
        Treatment = factor(rep(names(lst), each = length(dar_CORTCold)))) %>%
        mutate(Treatment = factor(Treatment, 
          levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
          labels = c("CORT-Cold" = paste0("CORT-Cold (n = ", n_list$delicata_Red_CORT_Cold, ")"),
                    "Control-Cold" = paste0("Control-Cold (n = ", n_list$delicata_Red_Control_Cold, ")"),
                    "CORT-Hot" = paste0("CORT-Hot (n = ", n_list$delicata_Red_CORT_Hot, ")"),
                    "Control-Hot" = paste0("Control-Hot (n = ", n_list$delicata_Red_Control_Hot, ")"))
          )
        ) %>%
      data.frame()
  } else if (data == "deli_blue"){
      if(var == "choice"){
        lst <- list('CORT-Cold' = probright_dbCORTCold,
                  'Control-Cold' = probright_dbControlCold,
                  'CORT-Hot' = probright_dbCORTHot,
                  'Control-Hot' = probright_dbControlHot)
      } else if(var == "slope"){
        lst <- list('CORT-Cold' = dab_CORTCold,
                  'Control-Cold' = dab_ControlCold,
                  'CORT-Hot' = dab_CORTHot,
                  'Control-Hot' = dab_ControlHot)
      }
      df_plotBD <- data.frame(Value = unlist(lst),
        Treatment = factor(rep(names(lst), each = length(dab_CORTCold)))) %>%
        mutate(Treatment = factor(Treatment, 
          levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
          labels = c("CORT-Cold" = paste0("CORT-Cold (n = ", n_list$delicata_Blue_CORT_Cold, ")"),
                    "Control-Cold" = paste0("Control-Cold (n = ", n_list$delicata_Blue_Control_Cold, ")"),
                    "CORT-Hot" = paste0("CORT-Hot (n = ", n_list$delicata_Blue_CORT_Hot, ")"),
                    "Control-Hot" = paste0("Control-Hot (n = ", n_list$delicata_Blue_Control_Hot, ")"))
          )
        ) %>%
      data.frame()    
  } else if(data == "guich_red"){
        if(var == "choice"){
          lst <- list('CORT-Cold' = probright_grCORTCold,
                  'Control-Cold' = probright_grControlCold,
                  'CORT-Hot' = probright_grCORTHot,
                  'Control-Hot' = probright_grControlHot)
        } else if(var == "slope"){
          lst <- list('CORT-Cold' = gar_CORTCold,
                    'Control-Cold' = gar_ControlCold,
                    'CORT-Hot' = gar_CORTHot,
                    'Control-Hot' = gar_ControlHot)
        }
        df_plotBD <- data.frame(Value = unlist(lst),
          Treatment = factor(rep(names(lst), each = length(gar_CORTCold)))) %>%
          mutate(Treatment = factor(Treatment, 
            levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
            labels = c("CORT-Cold" = paste0("CORT-Cold (n = ", n_list$guichenoti_Red_CORT_Cold, ")"),
                      "Control-Cold" = paste0("Control-Cold (n = ", n_list$guichenoti_Red_Control_Cold, ")"),
                      "CORT-Hot" = paste0("CORT-Hot (n = ", n_list$guichenoti_Red_CORT_Hot, ")"),
                      "Control-Hot" = paste0("Control-Hot (n = ", n_list$guichenoti_Red_Control_Hot, ")"))
            )
          ) %>%
      data.frame()
  } else if(data == "guich_blue"){
        if(var == "choice"){
          lst <- list('CORT-Cold' = probright_gbCORTCold,
                  'Control-Cold' = probright_gbControlCold,
                  'CORT-Hot' = probright_gbCORTHot,
                  'Control-Hot' = probright_gbControlHot)
        } else if(var == "slope"){
          lst <- list('CORT-Cold' = gab_CORTCold,
                    'Control-Cold' = gab_ControlCold,
                    'CORT-Hot' = gab_CORTHot,
                    'Control-Hot' = gab_ControlHot)
        }
        df_plotBD <- data.frame(Value = unlist(lst),
          Treatment = factor(rep(names(lst), each = length(gab_CORTCold)))) %>%
          mutate(Treatment = factor(Treatment, 
            levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot"),
            labels = c("CORT-Cold" = paste0("CORT-Cold (n = ", n_list$guichenoti_Blue_CORT_Cold, ")"),
                      "Control-Cold" = paste0("Control-Cold (n = ", n_list$guichenoti_Blue_Control_Cold, ")"),
                      "CORT-Hot" = paste0("CORT-Hot (n = ", n_list$guichenoti_Blue_CORT_Hot, ")"),
                      "Control-Hot" = paste0("Control-Hot (n = ", n_list$guichenoti_Blue_Control_Hot, ")"))
            )
          ) %>%
      data.frame()
  } else {
        stop("Data not valid")
  }
  return(df_plotBD)
}
####################
####################
# Function to create each df for the slopes figure (part of plots B, D fig-deli and fig-guich)
#' @title df_plotpoints
#' @param data to select the df, levels are: "deli_red", "deli_blue", "guich_red", "guich_blue"
#' @param var to select whether we use the intercept or the slope (choice/slope)
df_plotpoints <- function(data, var){
  if (data == "deli_red"){
    if(var == "choice"){
      prev_BD2 <- fig_dar_choiceviol_df
    } else if(var == "slope"){
      prev_BD2 <- fig_dar_sloviol_df
    }
  } else if (data == "deli_blue"){
    if(var == "choice"){
      prev_BD2 <- fig_dab_choiceviol_df
    } else if(var == "slope"){
      prev_BD2 <- fig_dab_sloviol_df
    }
  } else if (data == "guich_red"){
    if(var == "choice"){
      prev_BD2 <- fig_gar_choiceviol_df
    } else if(var == "slope"){
      prev_BD2 <- fig_gar_sloviol_df
    }
  }else if (data == "guich_blue"){
    if(var == "choice"){
      prev_BD2 <- fig_gab_choiceviol_df
    } else if(var == "slope"){
      prev_BD2 <- fig_gab_sloviol_df
    }
  } else {
        stop("Data not valid")
  }
df_BD2 <- prev_BD2 %>%
  group_by(Treatment) %>%
  summarize(
    Mean = mean(Value),
    SE = (sd(Value)/sqrt(length(Value))),
    SD = sd(Value)
  ) %>%
  ungroup() %>%
  data.frame()
  return(df_BD2)
}
####################
####################
# Function to create the plot per species (A-B or C-D) for fig-results
#' @title plotting
#' @param sp to select the species for the labels ("deli"/"guich")
#' @param col o select the colour ("red"/"blue")
#' @param df_violin_1 to select the df for the violin plot for choice (plot A, D)
#' @param df_points_1 to select the df for the points and geom_bars for choice (plot A, D)
#' @param df_violin_2 to select the df for the violin plot for slopes (plot B, E)
#' @param df_points_2 to select the df for the points and geom_bars for slopes (plot B, E)
#' @param df_prob to select the df for probability (plot C, F)
plotting <- function(sp, col, df_violin_1, df_points_1, df_violin_2, df_points_2, df_prob){
  # Specify labels depending on species and relevel the factor treatment for the legend
  if(sp == "deli"){
    if (col == "red"){
      lab <- c("A", "B", "C")
      custom_values <- c("CORT-Cold (n = 5)" = "#00008B",
                        "Control-Cold (n = 6)" = "#68bde1", 
                        "CORT-Hot (n = 5)" = "#b50101", 
                        "Control-Hot (n = 5)" = "#fa927d")
      custom_breaks <- c("Control-Hot (n = 5)",
                        "CORT-Hot (n = 5)",
                        "Control-Cold (n = 6)",
                        "CORT-Cold (n = 5)"
      )
      img <- readPNG("./Others/Red.png")

    } else if (col == "blue"){
      lab <- c("D", "E", "F")
            custom_values <- c("CORT-Cold (n = 6)" = "#00008B",
                        "Control-Cold (n = 6)" = "#68bde1", 
                        "CORT-Hot (n = 6)" = "#b50101", 
                        "Control-Hot (n = 5)" = "#fa927d")
      custom_breaks <- c("Control-Hot (n = 5)",
                        "CORT-Hot (n = 6)",
                        "Control-Cold (n = 6)",
                        "CORT-Cold (n = 6)"
      )
      img <- readPNG("./Others/Blue.png")
    } else {
      stop("Colour not valid")
    }
  } else if (sp == "guich"){
      if (col == "red"){
      lab <- c("A", "B", "C")
      custom_values <- c("CORT-Cold (n = 5)" = "#00008B",
                        "Control-Cold (n = 4)" = "#68bde1", 
                        "CORT-Hot (n = 5)" = "#b50101", 
                        "Control-Hot (n = 5)" = "#fa927d")
      custom_breaks <- c("Control-Hot (n = 5)",
                        "CORT-Hot (n = 5)",
                        "Control-Cold (n = 4)",
                        "CORT-Cold (n = 5)"
      )
      img <- readPNG("./Others/Red.png")
    } else if (col == "blue"){
      lab <- c("D", "E", "F")
      custom_values <- c("CORT-Cold (n = 5)" = "#00008B",
                        "Control-Cold (n = 3)" = "#68bde1", 
                        "CORT-Hot (n = 5)" = "#b50101", 
                        "Control-Hot (n = 5)" = "#fa927d")
      custom_breaks <- c("Control-Hot (n = 5)",
                        "CORT-Hot (n = 5)",
                        "Control-Cold (n = 3)",
                        "CORT-Cold (n = 5)"
      )
      img <- readPNG("./Others/Blue.png")
    } else {
      stop("Colour not valid")
    }
  } else {
    stop("Species not valid")
  }
  # First part of the plot (A, C), the probabilities of choosing right the first trial
  plot1 <- ggplot(df_violin_1, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_flat_violin(alpha = 0.5) +
  scale_fill_manual(values = custom_values,
                    breaks = custom_breaks) +
  geom_point(data = df_points_1, aes(y = Mean, x = Treatment), position = position_dodge(width = 0.75), color = "black", fill = "black", size = 3) +
  geom_segment(data = df_points_1, aes(y = Mean - SD, yend = Mean + SD, x = Treatment, xend = Treatment), size = 1.5, color = "black") +
  geom_hline(yintercept = 0.33, linetype = "dashed", color = "black") +
  ylim(min(df_violin_1$Value), max(df_violin_1$Value)) +
  coord_flip() +
  theme_classic() +
  labs(y = "Decision first trial", x = "Density") +
  theme(
    plot.margin = margin(7, 5.5, 5.5, 7, unit = "mm"),
    axis.title = element_text(size = 11, family = "Times"),
    axis.text.x = element_text(size = 9, family = "Times"),
    axis.text.y = element_blank(),
    legend.position = "none"
  ) 
  # Second part of the plot (B, E), the estimated slopes per treatment
  plot2 <- ggplot(df_violin_2, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_flat_violin(alpha = 0.5) +
  scale_fill_manual(values = custom_values,
                    breaks = custom_breaks) +
  geom_point(data = df_points_2, aes(y = Mean, x = Treatment), position = position_dodge(width = 0.75), color = "black", fill = "black", size = 3) +
  geom_segment(data = df_points_2, aes(y = Mean - SD, yend = Mean + SD, x = Treatment, xend = Treatment), size = 1.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ylim(-0.1, max(df_violin_2$Value)) +
  coord_flip() +
  theme_classic() +
  labs(y = "Slope estimates", x = "Density") +
  theme(
    plot.margin = margin(7, 5.5, 5.5, 7, unit = "mm"),
    axis.title = element_text(size = 11, family = "Times"),
    axis.text.x = element_text(size = 9, family = "Times"),
    axis.text.y = element_blank(),
    legend.position = "none"
  )
  # Third part of the plot (C, F), the estimated probability per treatments over trial
  plot3 <- ggplot(df_prob, aes(x = Trial, y = Mean_Predicted_prob, color = Treatment)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = custom_values,
                    breaks = custom_breaks) +
  geom_ribbon(aes(ymin = Mean_Predicted_prob - SE_Predicted_prob, ymax = Mean_Predicted_prob + SE_Predicted_prob, fill = Treatment), color = NA, alpha = 0.075) + 
  scale_fill_manual(values = custom_values,
                    breaks = custom_breaks) +
  theme_classic() +
  labs(y = "Probability of correct choice", x = "Trial") +
  theme(
    plot.margin = margin(7, 5.5, 5.5, 7, unit = "mm"),
    axis.title = element_text(size = 11, family = "Times"),
    axis.text = element_text(size = 9, family = "Times"),
    legend.position = "right",
    legend.title = element_text(size = 10, family = "Times"),
    legend.text = element_text(size = 9, family = "Times")
    )
  #
  #
  # Combine them
  plot <- plot_grid(plot1, plot2, plot3, labels = lab, nrow = 1, rel_widths = c(0.23, 0.23, 0.54), label_x = 0, label_y = 0.95) +
    annotation_custom(rasterGrob(img), xmin = 0.8, xmax = 0.98, ymin = 0.8, ymax = 0.98)
  return(plot)
}