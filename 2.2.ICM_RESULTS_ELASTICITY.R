### Integrated Community Model (ICM)
### Elk elasticity analysis
### Last updated: July 26, 2026
### xprockox@gmail.com

################################################################################
#########------------------- Load packages ---------------------################
################################################################################

library(tidyverse)
library(MCMCvis)
library(coda)
library(popbio)
library(cowplot)
library(stringr)

################################################################################
############------------------ Load results --------------------################
################################################################################

# load model results
load("data/outputs/ICM_parallel_output_2026-07-24.RData")

################################################################################
############--------------- Load helper functions --------------################
################################################################################

source("ICM_RESULTS_HELPER_FUNCTIONS.R")

################################################################################
#########--------------- Clean posterior samples ---------------################
################################################################################

cleaned <- clean_mcmc_samples(icm_samples)

icm_clean <- cleaned$icm_clean
post_mat <- cleaned$post_mat
bad_summary <- cleaned$bad_summary

# Review posterior columns removed because they contained nonfinite values
bad_summary

################################################################################
#########------------------ Vital-rate labels -------------------###############
################################################################################

elk_vrate_map <- c(
  elk_s_c = "Calf survival",
  elk_s_ya = "Young adult survival",
  elk_s_oa = "Old adult survival",
  elk_p_13 = "Young-to-old transition",
  elk_f_ya = "Young adult fecundity",
  elk_f_oa = "Old adult fecundity"
)

################################################################################
#########---------------- Vital-rate summaries -----------------################
################################################################################

elk_vrates <- extract_indexed_summary(
  mcmc_obj = icm_clean,
  params = c(
    "elk_s_c",
    "elk_s_ya",
    "elk_s_oa",
    "elk_p_13",
    "elk_f_ya",
    "elk_f_oa"
  ),
  years = community_years,
  label_map = elk_vrate_map,
  fecundity_params = c(
    "elk_f_ya",
    "elk_f_oa"
  )
)

################################################################################
#########---------------- Projection matrices ------------------################
################################################################################

elk_rates_wide <- elk_vrates %>%
  select(
    year,
    parameter_label,
    mean
  ) %>%
  pivot_wider(
    names_from = parameter_label,
    values_from = mean
  ) %>%
  drop_na(
    `Calf survival`,
    `Young adult survival`,
    `Old adult survival`,
    `Young-to-old transition`,
    `Young adult fecundity`,
    `Old adult fecundity`
  ) %>%
  arrange(year)

################################################################################
#########--------------- Annual elasticity analysis ------------################
################################################################################

elk_elasticity_df <- bind_rows(
  lapply(
    seq_len(nrow(elk_rates_wide)),
    function(i) {
      
      s_c <- elk_rates_wide$`Calf survival`[i]
      
      s_ya <- elk_rates_wide$`Young adult survival`[i]
      
      s_oa <- elk_rates_wide$`Old adult survival`[i]
      
      p_13 <- elk_rates_wide$`Young-to-old transition`[i]
      
      f_ya <- elk_rates_wide$`Young adult fecundity`[i]
      
      f_oa <- elk_rates_wide$`Old adult fecundity`[i]
      
      # Female elk projection matrix
      A <- matrix(
        c(
          0,
          f_ya * s_c,
          f_oa * s_c,
          
          s_ya,
          s_ya * (1 - p_13),
          0,
          
          0,
          s_ya * p_13,
          s_oa
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(
          to_stage = c(
            "Yearling",
            "Young adult",
            "Old adult"
          ),
          from_stage = c(
            "Yearling",
            "Young adult",
            "Old adult"
          )
        )
      )
      
      elasticity_matrix <- popbio::elasticity(A)
      
      tibble(
        year = elk_rates_wide$year[i],
        
        lambda = popbio::lambda(A),
        
        yearling_from_young_adult =
          elasticity_matrix[1, 2],
        
        yearling_from_old_adult =
          elasticity_matrix[1, 3],
        
        yearling_to_young_adult =
          elasticity_matrix[2, 1],
        
        young_adult_retention =
          elasticity_matrix[2, 2],
        
        young_adult_to_old_adult =
          elasticity_matrix[3, 2],
        
        old_adult_retention =
          elasticity_matrix[3, 3]
      )
    }
  )
)

################################################################################
#########---------------- Long-format elasticities -------------################
################################################################################

elk_elasticity_long <- elk_elasticity_df %>%
  pivot_longer(
    cols = c(
      yearling_from_young_adult,
      yearling_from_old_adult,
      yearling_to_young_adult,
      young_adult_retention,
      young_adult_to_old_adult,
      old_adult_retention
    ),
    names_to = "transition",
    values_to = "elasticity"
  ) %>%
  mutate(
    transition = factor(
      transition,
      levels = c(
        "yearling_from_young_adult",
        "yearling_from_old_adult",
        "yearling_to_young_adult",
        "young_adult_retention",
        "young_adult_to_old_adult",
        "old_adult_retention"
      ),
      labels = c(
        "Yearling recruitment from young adults",
        "Yearling recruitment from old adults",
        "Yearling to young adult",
        "Young adult survival without transition",
        "Young adult to old adult",
        "Old adult survival"
      )
    )
  )

################################################################################
#########---------------- Elasticity summaries -----------------################
################################################################################

elk_mean_elasticity <- elk_elasticity_long %>%
  group_by(transition) %>%
  summarise(
    mean_elasticity = mean(
      elasticity,
      na.rm = TRUE
    ),
    
    median_elasticity = median(
      elasticity,
      na.rm = TRUE
    ),
    
    minimum_elasticity = min(
      elasticity,
      na.rm = TRUE
    ),
    
    maximum_elasticity = max(
      elasticity,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    desc(mean_elasticity)
  )

elk_lambda_summary <- elk_elasticity_df %>%
  summarise(
    mean_lambda = mean(
      lambda,
      na.rm = TRUE
    ),
    
    median_lambda = median(
      lambda,
      na.rm = TRUE
    ),
    
    minimum_lambda = min(
      lambda,
      na.rm = TRUE
    ),
    
    maximum_lambda = max(
      lambda,
      na.rm = TRUE
    ),
    
    proportion_lambda_above_one = mean(
      lambda > 1,
      na.rm = TRUE
    )
  )

################################################################################
#########----------------------- Plots --------------------------################
################################################################################

lambda_plot <- ggplot(
  elk_elasticity_df,
  aes(
    x = year,
    y = lambda
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = 2
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 2
  ) +
  theme_classic() +
  labs(
    x = "Year",
    y = expression(lambda),
    title = "Annual elk population growth rate",
    subtitle = paste(
      "Calculated from annual posterior mean vital rates;",
      "the dashed line indicates lambda = 1"
    )
  )


elasticity_bar_plot <- elk_mean_elasticity %>%
  ggplot(
    aes(
      x = reorder(
        transition,
        mean_elasticity
      ),
      y = mean_elasticity
    )
  ) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean elasticity",
    title = "Mean elasticity of elk population growth",
    subtitle = "Elasticity to projection-matrix transitions"
  )


elasticity_time_plot <- ggplot(
  elk_elasticity_long,
  aes(
    x = year,
    y = elasticity,
    group = transition
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  facet_wrap(
    ~transition,
    scales = "free_y"
  ) +
  theme_classic() +
  labs(
    x = "Year",
    y = "Elasticity",
    title = "Annual elk elasticities by matrix transition"
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 8
    )
  )


elasticity_combo <- plot_grid(
  lambda_plot,
  elasticity_bar_plot,
  ncol = 1,
  rel_heights = c(
    1,
    1.1
  ),
  labels = c(
    "A",
    "B"
  )
)

################################################################################
#########--------------------- View results ---------------------###############
################################################################################

# Annual projection-matrix results
elk_elasticity_df

# Mean elasticity ranking
elk_mean_elasticity

# Summary of annual lambda values
elk_lambda_summary

# Plots
lambda_plot
elasticity_bar_plot
elasticity_time_plot
elasticity_combo

################################################################################
################################################################################
################################################################################