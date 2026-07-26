### Integrated Community Model (ICM)
### Results helper functions
### Last updated: July 26, 2026
### xprockox@gmail.com

################################################################################
#########------------------ Helper functions --------------------###############
################################################################################

# clean MCMC samples
clean_mcmc_samples <- function(chain_samples) {
  
  # Accept either a list of mcmc chains or an existing mcmc.list
  if (inherits(chain_samples, "mcmc.list")) {
    chain_list <- unclass(chain_samples)
  } else {
    chain_list <- chain_samples
  }
  
  post_mat_raw <- as.matrix(mcmc.list(chain_list))
  
  finite_cols <- apply(
    post_mat_raw,
    2,
    function(x) all(is.finite(x))
  )
  
  bad_cols <- which(!finite_cols)
  
  bad_summary <- tibble(
    parameter = colnames(post_mat_raw)[bad_cols],
    n_bad = vapply(
      bad_cols,
      function(i) sum(!is.finite(post_mat_raw[, i])),
      numeric(1)
    ),
    n_total = nrow(post_mat_raw)
  )
  
  good_cols <- colnames(post_mat_raw)[finite_cols]
  
  if (length(good_cols) == 0) {
    stop("No completely finite posterior columns were found.")
  }
  
  icm_clean <- mcmc.list(
    lapply(chain_list, function(ch) {
      mcmc(as.matrix(ch)[, good_cols, drop = FALSE])
    })
  )
  
  list(
    icm_clean = icm_clean,
    post_mat = as.matrix(icm_clean),
    bad_summary = bad_summary
  )
}

################################################################################
#########-------------- Validation helper functions -------------###############
################################################################################

# Extract posterior summaries for indexed abundance parameters
extract_abundance_summary <- function(mcmc_obj,
                                      params,
                                      years,
                                      stage_map) {
  
  # Keep only parameters that are present in the posterior samples
  available_params <- params[
    vapply(
      params,
      function(x) {
        any(str_detect(
          colnames(as.matrix(mcmc_obj)),
          paste0("^", x, "\\[")
        ))
      },
      logical(1)
    )
  ]
  
  missing_params <- setdiff(params, available_params)
  
  if (length(missing_params) > 0) {
    message(
      "The following abundance parameters were not found and will be skipped:\n",
      paste(missing_params, collapse = ", ")
    )
  }
  
  if (length(available_params) == 0) {
    return(NULL)
  }
  
  MCMCsummary(
    mcmc_obj,
    params = available_params
  ) %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    mutate(
      parameter_base = str_extract(parameter, "^[^\\[]+"),
      year_index = as.integer(
        str_extract(parameter, "(?<=\\[)\\d+(?=\\])")
      )
    ) %>%
    transmute(
      stage = recode(parameter_base, !!!stage_map),
      year = years[year_index],
      mean = mean,
      low = `2.5%`,
      high = `97.5%`
    )
}


# Convert observed abundance data to long format
make_observed_long <- function(data,
                               year_col,
                               cols_map) {
  
  missing_cols <- setdiff(
    c(year_col, names(cols_map)),
    names(data)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "The following columns are absent from the observed data:\n",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  data %>%
    select(all_of(c(year_col, names(cols_map)))) %>%
    rename(year = all_of(year_col)) %>%
    pivot_longer(
      cols = -year,
      names_to = "stage",
      values_to = "observed"
    ) %>%
    mutate(
      stage = recode(stage, !!!cols_map)
    )
}


# Plot posterior abundance estimates against observed values
plot_validation <- function(summary_df,
                            observed_df,
                            title,
                            y_label = "Abundance",
                            observed_stages = NULL) {
  
  if (is.null(summary_df) || nrow(summary_df) == 0) {
    return(NULL)
  }
  
  if (!is.null(observed_stages)) {
    observed_df <- observed_df %>%
      filter(stage %in% observed_stages)
  }
  
  ggplot(
    summary_df,
    aes(
      x = year,
      y = mean,
      group = stage
    )
  ) +
    geom_ribbon(
      aes(
        ymin = low,
        ymax = high,
        fill = stage
      ),
      alpha = 0.20
    ) +
    geom_line(
      linewidth = 0.9
    ) +
    geom_line(
      data = observed_df,
      aes(
        x = year,
        y = observed,
        group = stage
      ),
      inherit.aes = FALSE,
      color = "red",
      linetype = 2,
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    geom_point(
      data = observed_df,
      aes(
        x = year,
        y = observed
      ),
      inherit.aes = FALSE,
      color = "red",
      size = 2,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~stage,
      scales = "free_y"
    ) +
    theme_bw() +
    labs(
      x = "Year",
      y = y_label,
      title = title,
      subtitle = paste(
        "Posterior mean and 95% credible interval;",
        "observed values shown in red"
      )
    ) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold")
    )
}

################################################################################
#########--------------- Density helper functions ---------------###############
################################################################################

resolve_coef_keep <- function(coef_names, label_map, keep = "all") {
  
  pretty_names <- unname(label_map[coef_names])
  pretty_names <- pretty_names[!is.na(pretty_names)]
  
  if (length(keep) == 1) {
    
    if (identical(keep, "all")) {
      return(pretty_names)
    }
    
    if (identical(keep, "intercepts")) {
      return(
        pretty_names[
          str_detect(
            pretty_names,
            regex("intercept", ignore_case = TRUE)
          )
        ]
      )
    }
    
    if (identical(keep, "non_intercepts")) {
      return(
        pretty_names[
          !str_detect(
            pretty_names,
            regex("intercept", ignore_case = TRUE)
          )
        ]
      )
    }
    
    if (str_detect(keep, "^regex:")) {
      pattern <- str_remove(keep, "^regex:")
      
      return(
        pretty_names[
          str_detect(
            pretty_names,
            regex(pattern, ignore_case = TRUE)
          )
        ]
      )
    }
  }
  
  pretty_names[pretty_names %in% keep]
}


get_stage_coef_names <- function(model_specs,
                                 stage_name,
                                 post_mat,
                                 include_intercept = TRUE) {
  
  if (!stage_name %in% names(model_specs)) {
    stop("Stage not found in model_specs: ", stage_name)
  }
  
  stage_spec <- model_specs[[stage_name]]
  
  coef_names <- unname(stage_spec$terms)
  
  if (include_intercept) {
    coef_names <- c(stage_spec$intercept, coef_names)
  }
  
  missing_coefs <- setdiff(coef_names, colnames(post_mat))
  
  if (length(missing_coefs) > 0) {
    message(
      "The following coefficients are listed in the model specification ",
      "but are absent from post_mat and will be skipped:\n",
      paste(missing_coefs, collapse = ", ")
    )
  }
  
  intersect(coef_names, colnames(post_mat))
}


get_all_model_coef_names <- function(model_specs,
                                     post_mat,
                                     include_intercepts = TRUE) {
  
  coef_names <- unlist(
    lapply(model_specs, function(stage_spec) {
      
      out <- unname(stage_spec$terms)
      
      if (include_intercepts) {
        out <- c(stage_spec$intercept, out)
      }
      
      out
    }),
    use.names = FALSE
  )
  
  coef_names <- unique(coef_names)
  
  missing_coefs <- setdiff(coef_names, colnames(post_mat))
  
  if (length(missing_coefs) > 0) {
    message(
      "The following model coefficients are absent from post_mat ",
      "and will be skipped:\n",
      paste(missing_coefs, collapse = ", ")
    )
  }
  
  intersect(coef_names, colnames(post_mat))
}


make_coef_density_plot_stacked <- function(
    post_mat,
    coef_names,
    label_map,
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
    model_order = c(
      "Calf survival",
      "YA survival",
      "OA survival",
      "Wolf pup survival",
      "Wolf adult survival"
    ),
    scale_factor = 0.35) {
  
  coef_names <- intersect(coef_names, colnames(post_mat))
  
  if (length(coef_names) == 0) {
    stop("None of the requested coefficients were found in post_mat.")
  }
  
  missing_labels <- setdiff(coef_names, names(label_map))
  
  if (length(missing_labels) > 0) {
    stop(
      "The following coefficients do not have entries in label_map:\n",
      paste(missing_labels, collapse = ", ")
    )
  }
  
  keep_resolved <- resolve_coef_keep(
    coef_names = coef_names,
    label_map = label_map,
    keep = keep
  )
  
  if (length(keep_resolved) == 0) {
    stop(
      "The coefficient selector did not retain any coefficients. ",
      "Check the value supplied to 'keep'."
    )
  }
  
  df <- as_tibble(
    post_mat[, coef_names, drop = FALSE]
  ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
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
        str_detect(
          pretty,
          regex("intercept", ignore_case = TRUE)
        ) ~ "Intercept",
        
        str_detect(
          pretty,
          regex("grizzly effect", ignore_case = TRUE)
        ) ~ "Grizzly abundance",
        
        str_detect(
          pretty,
          regex("bison effect", ignore_case = TRUE)
        ) ~ "Bison abundance",
        
        str_detect(
          pretty,
          regex("cougar effect", ignore_case = TRUE)
        ) ~ "Cougar abundance",
        
        str_detect(
          pretty,
          regex("annual NPP effect", ignore_case = TRUE)
        ) ~ "Annual NPP",
        
        str_detect(
          pretty,
          regex("browndown effect", ignore_case = TRUE)
        ) ~ "Browndown",
        
        str_detect(
          pretty,
          regex("PDSI effect", ignore_case = TRUE)
        ) ~ "PDSI",
        
        str_detect(
          pretty,
          regex("winter severity effect", ignore_case = TRUE)
        ) ~ "Winter severity",
        
        str_detect(
          pretty,
          regex("density dependence effect", ignore_case = TRUE)
        ) ~ "Density dependence",
        
        str_detect(
          pretty,
          regex("wolf effect", ignore_case = TRUE)
        ) ~ "Wolf abundance",
        
        str_detect(
          pretty,
          regex("elk effect", ignore_case = TRUE)
        ) ~ "Elk abundance",
        
        TRUE ~ pretty
      ),
      
      effect_type = factor(
        effect_type,
        levels = effect_order
      ),
      
      model_group = factor(
        model_group,
        levels = model_order
      )
    ) %>%
    drop_na(model_group, effect_type)
  
  if (nrow(df) == 0) {
    stop(
      "No coefficients remained after assigning model groups and effect types."
    )
  }
  
  dens_df <- df %>%
    group_by(model_group, effect_type) %>%
    group_modify(~ {
      
      values <- .x$value[is.finite(.x$value)]
      
      if (length(unique(values)) < 2) {
        return(tibble())
      }
      
      d <- density(values)
      
      tibble(
        x = d$x,
        density = d$y
      )
    }) %>%
    ungroup() %>%
    mutate(
      effect_num = as.numeric(effect_type),
      y = effect_num + density * scale_factor
    )
  
  if (nrow(dens_df) == 0) {
    stop("Density estimation failed for all requested coefficients.")
  }
  
  label_df <- df %>%
    distinct(model_group, effect_type) %>%
    mutate(
      effect_num = as.numeric(effect_type)
    )
  
  ggplot() +
    geom_ribbon(
      data = dens_df,
      aes(
        x = x,
        ymin = effect_num,
        ymax = y,
        group = interaction(model_group, effect_type)
      ),
      fill = fill_color,
      alpha = 0.45
    ) +
    geom_line(
      data = dens_df,
      aes(
        x = x,
        y = y,
        group = interaction(model_group, effect_type)
      ),
      linewidth = 0.5
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    geom_hline(
      data = label_df,
      aes(yintercept = effect_num),
      color = "grey80",
      linewidth = 0.3
    ) +
    facet_wrap(
      ~model_group,
      nrow = 1,
      scales = "free_x"
    ) +
    scale_y_continuous(
      breaks = sort(unique(label_df$effect_num)),
      labels = levels(droplevels(label_df$effect_type)),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    theme_bw() +
    labs(
      x = "Posterior coefficient value",
      y = NULL,
      title = title
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}


make_single_coef_density_plot <- function(post_mat,
                                          parameter,
                                          fill_color,
                                          title,
                                          x_label = "Posterior coefficient value") {
  
  if (!parameter %in% colnames(post_mat)) {
    message(
      "Parameter '",
      parameter,
      "' was not found in post_mat. Plot was not created."
    )
    
    return(NULL)
  }
  
  tibble(
    value = post_mat[, parameter]
  ) %>%
    ggplot(aes(x = value)) +
    geom_density(
      fill = fill_color,
      alpha = 0.45
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    theme_bw() +
    labs(
      x = x_label,
      y = "Density",
      title = title
    )
}


################################################################################
#########------------------ Parameter utilities -----------------###############
################################################################################

# Find all indexed columns for a parameter and order them numerically
get_indexed_columns <- function(post_mat, parameter) {
  
  cols <- grep(
    paste0("^", parameter, "\\["),
    colnames(post_mat),
    value = TRUE
  )
  
  if (length(cols) == 0) {
    return(character(0))
  }
  
  indices <- as.integer(
    str_extract(cols, "(?<=\\[)\\d+(?=\\])")
  )
  
  cols[order(indices)]
}


# Extract posterior draws for an indexed parameter
get_indexed_draws <- function(post_mat, parameter) {
  
  cols <- get_indexed_columns(
    post_mat = post_mat,
    parameter = parameter
  )
  
  if (length(cols) == 0) {
    stop(
      "No posterior columns were found for parameter: ",
      parameter
    )
  }
  
  post_mat[, cols, drop = FALSE]
}


# Summarize an indexed posterior parameter by year
summarize_indexed_parameter <- function(post_mat,
                                        parameter,
                                        years,
                                        value_name) {
  
  draws <- get_indexed_draws(
    post_mat = post_mat,
    parameter = parameter
  )
  
  if (ncol(draws) != length(years)) {
    stop(
      "The number of posterior columns for ",
      parameter,
      " does not match the number of supplied years."
    )
  }
  
  tibble(
    year = years,
    !!value_name := colMeans(draws),
    low = apply(
      draws,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    high = apply(
      draws,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    )
  )
}


# Return TRUE when a scalar posterior coefficient exists
coefficient_available <- function(post_mat, coefficient) {
  coefficient %in% colnames(post_mat)
}


################################################################################
#########-------------- Build survival-abundance points ----------##############
################################################################################

make_survival_abundance_points <- function(post_mat,
                                           model_specs,
                                           predictor_name,
                                           predictor_specs,
                                           community_years,
                                           regression_years) {
  
  predictor_spec <- predictor_specs[[predictor_name]]
  
  if (is.null(predictor_spec)) {
    stop(
      "No predictor specification was found for: ",
      predictor_name
    )
  }
  
  abundance_summary <- summarize_indexed_parameter(
    post_mat = post_mat,
    parameter = predictor_spec$parameter,
    years = community_years,
    value_name = "abundance"
  ) %>%
    transmute(
      predictor_year = year,
      survival_year = year + 1,
      abundance,
      abundance_low = low,
      abundance_high = high
    )
  
  stage_points <- lapply(
    names(model_specs),
    function(stage_name) {
      
      stage_spec <- model_specs[[stage_name]]
      
      if (!predictor_name %in% names(stage_spec$terms)) {
        return(NULL)
      }
      
      coefficient <- unname(
        stage_spec$terms[[predictor_name]]
      )
      
      if (!coefficient_available(post_mat, coefficient)) {
        message(
          "Skipping ",
          stage_name,
          " × ",
          predictor_name,
          " because coefficient ",
          coefficient,
          " is absent from post_mat."
        )
        
        return(NULL)
      }
      
      survival_summary <- summarize_indexed_parameter(
        post_mat = post_mat,
        parameter = stage_spec$survival_parameter,
        years = community_years,
        value_name = "survival"
      ) %>%
        transmute(
          survival_year = year,
          survival,
          survival_low = low,
          survival_high = high
        )
      
      survival_summary %>%
        inner_join(
          abundance_summary,
          by = "survival_year"
        ) %>%
        filter(
          survival_year %in% regression_years
        ) %>%
        mutate(
          stage = stage_name,
          predictor = predictor_name
        )
    }
  )
  
  bind_rows(stage_points)
}

################################################################################
#########--------------- Build regression curves ----------------###############
################################################################################

make_survival_abundance_curve <- function(post_mat,
                                          model_specs,
                                          predictor_name,
                                          predictor_specs,
                                          point_df,
                                          length_out = 200) {
  
  predictor_spec <- predictor_specs[[predictor_name]]
  
  if (nrow(point_df) == 0) {
    return(tibble())
  }
  
  x_values <- seq(
    min(point_df$abundance, na.rm = TRUE),
    max(point_df$abundance, na.rm = TRUE),
    length.out = length_out
  )
  
  x_standardized <- (
    x_values - predictor_spec$mean
  ) / predictor_spec$sd
  
  stage_curves <- lapply(
    names(model_specs),
    function(stage_name) {
      
      stage_spec <- model_specs[[stage_name]]
      
      if (!predictor_name %in% names(stage_spec$terms)) {
        return(NULL)
      }
      
      coefficient <- unname(
        stage_spec$terms[[predictor_name]]
      )
      
      required_parameters <- c(
        stage_spec$intercept,
        coefficient
      )
      
      missing_parameters <- setdiff(
        required_parameters,
        colnames(post_mat)
      )
      
      if (length(missing_parameters) > 0) {
        return(NULL)
      }
      
      intercept_draws <- post_mat[
        ,
        stage_spec$intercept,
        drop = TRUE
      ]
      
      coefficient_draws <- post_mat[
        ,
        coefficient,
        drop = TRUE
      ]
      
      # All nonfocal standardized predictors are held at zero,
      # corresponding to their standardization means.
      eta <- outer(
        coefficient_draws,
        x_standardized
      ) + intercept_draws
      
      survival_draws <- plogis(eta)
      
      tibble(
        abundance = x_values,
        survival = colMeans(survival_draws),
        survival_low = apply(
          survival_draws,
          2,
          quantile,
          probs = 0.025,
          na.rm = TRUE
        ),
        survival_high = apply(
          survival_draws,
          2,
          quantile,
          probs = 0.975,
          na.rm = TRUE
        ),
        stage = stage_name,
        predictor = predictor_name
      )
    }
  )
  
  bind_rows(stage_curves)
}

################################################################################
#########----------------- Plotting function --------------------###############
################################################################################

plot_survival_abundance <- function(point_df,
                                    curve_df,
                                    x_label,
                                    title) {
  
  if (nrow(point_df) == 0 || nrow(curve_df) == 0) {
    message("No plot was created for: ", title)
    return(NULL)
  }
  
  ggplot(
    point_df,
    aes(
      x = abundance,
      y = survival
    )
  ) +
    geom_ribbon(
      data = curve_df,
      aes(
        x = abundance,
        ymin = survival_low,
        ymax = survival_high
      ),
      inherit.aes = FALSE,
      alpha = 0.20
    ) +
    geom_line(
      data = curve_df,
      aes(
        x = abundance,
        y = survival
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    geom_errorbar(
      aes(
        ymin = survival_low,
        ymax = survival_high
      ),
      width = 0,
      alpha = 0.55
    ) +
    geom_errorbarh(
      aes(
        xmin = abundance_low,
        xmax = abundance_high
      ),
      height = 0,
      alpha = 0.55
    ) +
    geom_point(
      shape = 21,
      fill = "gold",
      color = "black",
      size = 2.5
    ) +
    facet_wrap(
      ~stage,
      scales = "free_y"
    ) +
    theme_classic() +
    labs(
      x = paste0(x_label, " in year t - 1"),
      y = "Survival in year t",
      title = title,
      subtitle = paste(
        "Points show annual posterior abundance-survival estimates;",
        "curve shows the conditional regression effect with other predictors held at zero"
      )
    ) +
    theme(
      strip.text = element_text(face = "bold")
    )
}


make_conditional_points <- function(post_mat,
                                    model_specs,
                                    predictor_name,
                                    predictor_specs,
                                    abundance_df) {
  
  predictor_spec <- predictor_specs[[predictor_name]]
  
  bind_rows(lapply(names(model_specs), function(stage_name) {
    
    stage_spec <- model_specs[[stage_name]]
    
    if (!predictor_name %in% names(stage_spec$terms)) {
      return(NULL)
    }
    
    beta_name <- unname(stage_spec$terms[[predictor_name]])
    
    required <- c(stage_spec$intercept, beta_name)
    
    if (!all(required %in% colnames(post_mat))) {
      return(NULL)
    }
    
    intercept_draws <- post_mat[, stage_spec$intercept]
    beta_draws <- post_mat[, beta_name]
    
    bind_rows(lapply(seq_len(nrow(abundance_df)), function(i) {
      
      abundance_i <- abundance_df$abundance[i]
      
      abundance_std_i <- (
        abundance_i - predictor_spec$mean
      ) / predictor_spec$sd
      
      survival_draws <- plogis(
        intercept_draws +
          beta_draws * abundance_std_i
      )
      
      tibble(
        survival_year = abundance_df$survival_year[i],
        abundance = abundance_i,
        abundance_low = abundance_df$abundance_low[i],
        abundance_high = abundance_df$abundance_high[i],
        survival = mean(survival_draws),
        survival_low = quantile(
          survival_draws,
          0.025,
          na.rm = TRUE
        ),
        survival_high = quantile(
          survival_draws,
          0.975,
          na.rm = TRUE
        ),
        stage = stage_name,
        predictor = predictor_name
      )
    }))
  }))
}

# Build a dataframe of lagged annual abundance estimates
make_lagged_abundance_df <- function(post_mat,
                                     predictor_name,
                                     predictor_specs,
                                     community_years,
                                     regression_years) {
  
  predictor_spec <- predictor_specs[[predictor_name]]
  
  if (is.null(predictor_spec)) {
    stop(
      "No predictor specification was found for: ",
      predictor_name
    )
  }
  
  summarize_indexed_parameter(
    post_mat = post_mat,
    parameter = predictor_spec$parameter,
    years = community_years,
    value_name = "abundance"
  ) %>%
    transmute(
      predictor_year = year,
      survival_year = year + 1,
      abundance,
      abundance_low = low,
      abundance_high = high
    ) %>%
    filter(
      survival_year %in% regression_years
    )
}


################################################################################
#########--------------- Elasticity functions -------------------###############
################################################################################

# Extract annual posterior summaries for indexed parameters
extract_indexed_summary <- function(mcmc_obj,
                                    params,
                                    years,
                                    label_map = NULL,
                                    fecundity_params = NULL) {
  
  posterior_names <- colnames(
    as.matrix(mcmc_obj)
  )
  
  available_params <- params[
    vapply(
      params,
      function(parameter) {
        any(str_detect(
          posterior_names,
          paste0("^", parameter, "\\[")
        ))
      },
      logical(1)
    )
  ]
  
  missing_params <- setdiff(
    params,
    available_params
  )
  
  if (length(missing_params) > 0) {
    message(
      "The following parameters were not found and will be skipped:\n",
      paste(missing_params, collapse = ", ")
    )
  }
  
  if (length(available_params) == 0) {
    stop("None of the requested parameters were found.")
  }
  
  out <- MCMCsummary(
    mcmc_obj,
    params = available_params
  ) %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    rename(
      low = `2.5%`,
      high = `97.5%`
    ) %>%
    mutate(
      parameter_base = str_extract(
        parameter,
        "^[^\\[]+"
      ),
      
      year_index = as.integer(
        str_extract(
          parameter,
          "(?<=\\[)\\d+(?=\\])"
        )
      )
    )
  
  if (!is.null(label_map)) {
    out <- out %>%
      mutate(
        parameter_label = recode(
          parameter_base,
          !!!label_map
        )
      )
  } else {
    out <- out %>%
      mutate(
        parameter_label = parameter_base
      )
  }
  
  out <- out %>%
    mutate(
      year = case_when(
        !is.null(fecundity_params) &
          parameter_base %in% fecundity_params ~
          years[-length(years)][year_index],
        
        TRUE ~ years[year_index]
      )
    )
  
  out
}


################################################################################
#########------------ Qualitative summary functions -------------###############
################################################################################

# Classify posterior evidence according to nested credible intervals
classify_evidence <- function(draws) {
  
  draws <- draws[is.finite(draws)]
  
  if (length(draws) == 0) {
    stop("No finite posterior draws were supplied.")
  }
  
  ci95 <- quantile(
    draws,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  ci80 <- quantile(
    draws,
    probs = c(0.10, 0.90),
    na.rm = TRUE
  )
  
  ci50 <- quantile(
    draws,
    probs = c(0.25, 0.75),
    na.rm = TRUE
  )
  
  excludes_zero_95 <- ci95[1] > 0 | ci95[2] < 0
  excludes_zero_80 <- ci80[1] > 0 | ci80[2] < 0
  excludes_zero_50 <- ci50[1] > 0 | ci50[2] < 0
  
  posterior_median <- median(draws)
  
  tibble(
    mean = mean(draws),
    median = posterior_median,
    
    ci95_low = unname(ci95[1]),
    ci95_high = unname(ci95[2]),
    
    ci80_low = unname(ci80[1]),
    ci80_high = unname(ci80[2]),
    
    ci50_low = unname(ci50[1]),
    ci50_high = unname(ci50[2]),
    
    probability_positive = mean(draws > 0),
    probability_negative = mean(draws < 0),
    
    direction = case_when(
      posterior_median > 0 ~ "Positive",
      posterior_median < 0 ~ "Negative",
      TRUE ~ "Neutral"
    ),
    
    evidence = case_when(
      excludes_zero_95 ~ "Strong",
      !excludes_zero_95 & excludes_zero_80 ~ "Moderate",
      !excludes_zero_80 & excludes_zero_50 ~ "Weak",
      TRUE ~ "Little/none"
    )
  )
}


# Extract summaries for indexed abundance parameters
extract_abundance_summary <- function(mcmc_obj,
                                      params,
                                      years,
                                      stage_map) {
  
  posterior_names <- colnames(as.matrix(mcmc_obj))
  
  available_params <- params[
    vapply(
      params,
      function(parameter) {
        any(str_detect(
          posterior_names,
          paste0("^", parameter, "\\[")
        ))
      },
      logical(1)
    )
  ]
  
  missing_params <- setdiff(
    params,
    available_params
  )
  
  if (length(missing_params) > 0) {
    message(
      "The following abundance parameters were not found and will be skipped:\n",
      paste(missing_params, collapse = ", ")
    )
  }
  
  if (length(available_params) == 0) {
    stop("None of the requested abundance parameters were found.")
  }
  
  MCMCsummary(
    mcmc_obj,
    params = available_params
  ) %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    mutate(
      parameter_base = str_extract(
        parameter,
        "^[^\\[]+"
      ),
      year_index = as.integer(
        str_extract(
          parameter,
          "(?<=\\[)\\d+(?=\\])"
        )
      )
    ) %>%
    transmute(
      stage = recode(
        parameter_base,
        !!!stage_map
      ),
      year = years[year_index],
      mean = mean,
      low = `2.5%`,
      high = `97.5%`
    )
}


# Format evidence results for console output
format_evidence_lines <- function(evidence_df,
                                  interval = c("95", "80", "50")) {
  
  interval <- match.arg(interval)
  
  interval_low <- paste0("ci", interval, "_low")
  interval_high <- paste0("ci", interval, "_high")
  
  if (nrow(evidence_df) == 0) {
    return("None")
  }
  
  evidence_df %>%
    filter(
      !str_detect(
        label,
        regex("intercept", ignore_case = TRUE)
      )
    ) %>%
    mutate(
      text = paste0(
        label,
        ": median = ",
        round(median, 3),
        ", ",
        interval,
        "% CI [",
        round(.data[[interval_low]], 3),
        ", ",
        round(.data[[interval_high]], 3),
        "]",
        ", direction = ",
        direction
      )
    ) %>%
    pull(text) %>%
    paste0("- ", ., collapse = "\n")
}
