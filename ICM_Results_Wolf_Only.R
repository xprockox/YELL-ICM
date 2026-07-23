### Integrated Community Model (ICM)
### Results exploration
### Last updated: July 9, 2026
### xprockox@gmail.com

################################################################################
#########------------------- Load packages ---------------------################
################################################################################

library(tidyverse)
library(MCMCvis)
library(coda)
library(cowplot)
library(popbio)
library(stringr)

################################################################################
############------------------ Load data ---------------------##################
################################################################################

# load model results
load("data/outputs/ICM_environment_2026-04-17.RData")

resolve_coef_keep <- function(coef_names, label_map, keep = "all") {
  pretty_names <- unname(label_map[coef_names])
  pretty_names <- pretty_names[!is.na(pretty_names)]
  
  if (length(keep) == 1) {
    if (identical(keep, "all")) {
      return(pretty_names)
    }
    
    if (identical(keep, "intercepts")) {
      return(pretty_names[str_detect(pretty_names, regex("intercept", ignore_case = TRUE))])
    }
    
    if (identical(keep, "non_intercepts")) {
      return(pretty_names[!str_detect(pretty_names, regex("intercept", ignore_case = TRUE))])
    }
    
    if (str_detect(keep, "^regex:")) {
      pattern <- str_remove(keep, "^regex:")
      return(pretty_names[str_detect(pretty_names, regex(pattern, ignore_case = TRUE))])
    }
  }
  
  pretty_names[pretty_names %in% keep]
}

make_coef_density_plot_stacked <- function(post_mat, coef_names, label_map,
                                           keep = "all",
                                           fill_color,
                                           title,
                                           effect_order = c(
                                             "Intercept",
                                             "Wolf abundance",
                                             "Grizzly abundance",
                                             "Cougar abundance",
                                             "Elk abundance",
                                             "Bison abundance",
                                             "Density dependence",
                                             "Winter severity",
                                             "Annual NPP",
                                             "Browndown",
                                             "PDSI"
                                           ),
                                           model_order = c("Calf survival", "YA survival", "OA survival",
                                                           "Wolf pup survival", "Wolf adult survival"),
                                           scale_factor = 0.35) {
  
  keep_resolved <- resolve_coef_keep(
    coef_names = coef_names,
    label_map = label_map,
    keep = keep
  )
  
  df <- as_tibble(post_mat[, coef_names, drop = FALSE]) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      pretty = recode(parameter, !!!label_map)
    ) %>%
    filter(pretty %in% keep_resolved) %>%
    mutate(
      model_group = case_when(
        str_detect(pretty, "^Calf survival") ~ "Calf survival",
        str_detect(pretty, "^YA survival") ~ "YA survival",
        str_detect(pretty, "^OA survival") ~ "OA survival",
        str_detect(pretty, "^Wolf pup survival") ~ "Wolf pup survival",
        str_detect(pretty, "^Wolf adult survival") ~ "Wolf adult survival",
        TRUE ~ "Other"
      ),
      effect_type = case_when(
        str_detect(pretty, regex("intercept", ignore_case = TRUE)) ~ "Intercept",
        str_detect(pretty, regex("grizzly effect", ignore_case = TRUE)) ~ "Grizzly abundance",
        str_detect(pretty, regex("bison effect", ignore_case = TRUE)) ~ "Bison abundance",
        str_detect(pretty, regex("cougar effect", ignore_case = TRUE)) ~ "Cougar abundance",
        str_detect(pretty, regex("annual NPP effect", ignore_case = TRUE)) ~ "Annual NPP",
        str_detect(pretty, regex("browndown effect", ignore_case = TRUE)) ~ "Browndown",
        str_detect(pretty, regex("PDSI effect", ignore_case = TRUE)) ~ "PDSI",
        str_detect(pretty, regex("winter severity effect", ignore_case = TRUE)) ~ "Winter severity",
        str_detect(pretty, regex("density dependence effect", ignore_case = TRUE)) ~ "Density dependence",
        str_detect(pretty, regex("wolf effect", ignore_case = TRUE)) ~ "Wolf abundance",
        str_detect(pretty, regex("elk effect", ignore_case = TRUE)) ~ "Elk abundance",
        TRUE ~ pretty
      )
    ) %>%
    mutate(
      effect_type = factor(effect_type, levels = effect_order),
      model_group = factor(model_group, levels = model_order)
    ) %>%
    drop_na(model_group, effect_type)
  
  dens_df <- df %>%
    group_by(model_group, effect_type) %>%
    group_modify(~ {
      d <- density(.x$value, na.rm = TRUE)
      tibble(x = d$x, density = d$y)
    }) %>%
    ungroup() %>%
    mutate(
      effect_num = as.numeric(effect_type),
      y = effect_num + density * scale_factor
    )
  
  label_df <- df %>%
    distinct(model_group, effect_type) %>%
    mutate(effect_num = as.numeric(effect_type))
  
  ggplot() +
    geom_ribbon(
      data = dens_df,
      aes(x = x, ymin = effect_num, ymax = y, group = interaction(model_group, effect_type)),
      fill = fill_color,
      alpha = 0.45
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = y, group = interaction(model_group, effect_type)),
      linewidth = 0.5
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(
      data = label_df,
      aes(yintercept = effect_num),
      color = "grey80",
      linewidth = 0.3
    ) +
    facet_wrap(~model_group, nrow = 1, scales = "free_x") +
    scale_y_continuous(
      breaks = label_df$effect_num,
      labels = label_df$effect_type,
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    theme_bw() +
    labs(
      x = "Posterior value",
      y = NULL,
      title = title
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}


old_wolf_names <- c(
  "beta1_calfSurv_wolfN",
  "beta1_yaSurv_wolfN",
  "beta1_oaSurv_wolfN"
)

# Combine every chain
old_post_mat <- do.call(
  rbind,
  lapply(icm_mod$samples, as.matrix)
)

sapply(
  old_wolf_names,
  function(x) sum(is.finite(old_post_mat[, x]))
)

old_label_map <- c(
  beta1_calfSurv_wolfN = "Calf survival: wolf effect",
  beta1_yaSurv_wolfN   = "YA survival: wolf effect",
  beta1_oaSurv_wolfN   = "OA survival: wolf effect"
)

old_wolf_plot <- make_coef_density_plot_stacked(
  post_mat = old_post_mat,
  coef_names = old_wolf_names,
  label_map = old_label_map,
  keep = "all",
  fill_color = "#236192",
  title = "Posterior distributions of wolf effects on elk survival",
  effect_order = "Wolf abundance",
  model_order = c(
    "Calf survival",
    "YA survival",
    "OA survival"
  )
)

old_wolf_plot