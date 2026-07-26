### Aspen Data Management for ICM
### Last updated: June 11, 2026
### xprockox@gmail.com

#################################################################################################
# packages

library(tidyverse)
library(cowplot)

#################################################################################################
# settings

height_threshold <- 120

#################################################################################################
# data import

dat <- read.csv('data/2024-02-25_AspenData_Clean.csv')

#################################################################################################
# data management

# first calculate plot-level counts of stems above the threshold within each year
plot_over_threshold_summary <- dat %>%
  group_by(year, plot) %>%
  summarise(
    n_over_threshold = sum(height_cm > height_threshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(year) %>%
  summarise(
    mean_n_stems_over_threshold_per_plot = mean(n_over_threshold, na.rm = TRUE),
    median_n_stems_over_threshold_per_plot = median(n_over_threshold, na.rm = TRUE),
    .groups = "drop"
  )

annual_aspen_metrics <- dat %>%
  group_by(year) %>%
  summarise(
    total_stems = n(),
    n_plots = n_distinct(plot),
    mean_stems_per_plot = total_stems / n_plots,
    
    n_new = sum(new == 1, na.rm = TRUE),
    prop_new = mean(new == 1, na.rm = TRUE),
    
    mean_height = mean(height_cm, na.rm = TRUE),
    median_height = median(height_cm, na.rm = TRUE),
    q90_height = quantile(height_cm, 0.90, na.rm = TRUE),
    
    n_stems_over_threshold = sum(height_cm > height_threshold, na.rm = TRUE),
    prop_over_threshold = mean(height_cm > height_threshold, na.rm = TRUE),
    
    mean_cag = mean(cag_cm, na.rm = TRUE),
    median_cag = median(cag_cm, na.rm = TRUE),
    
    prop_browsed = mean(browse_winter == 1, na.rm = TRUE),
    
    n_over_threshold_unbrowsed = sum(
      height_cm > height_threshold & browse_winter == 0,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  left_join(plot_over_threshold_summary, by = "year")

annual_aspen_plot_metrics <- dat %>%
  group_by(year, plot) %>%
  summarise(
    n_stems = n(),
    n_new = sum(new == 1, na.rm = TRUE),
    n_over_threshold = sum(height_cm > height_threshold, na.rm = TRUE),
    prop_over_threshold = mean(height_cm > height_threshold, na.rm = TRUE),
    median_height = median(height_cm, na.rm = TRUE),
    prop_browsed = mean(browse_winter == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(year) %>%
  summarise(
    mean_plot_stems = mean(n_stems, na.rm = TRUE),
    median_plot_stems = median(n_stems, na.rm = TRUE),
    mean_plot_over_threshold = mean(n_over_threshold, na.rm = TRUE),
    median_plot_over_threshold = median(n_over_threshold, na.rm = TRUE),
    prop_plots_with_over_threshold = mean(n_over_threshold > 0, na.rm = TRUE),
    prop_plots_with_5plus_over_threshold = mean(n_over_threshold >= 5, na.rm = TRUE),
    mean_plot_median_height = mean(median_height, na.rm = TRUE),
    mean_plot_browse = mean(prop_browsed, na.rm = TRUE),
    .groups = "drop"
  )

#################################################################################################
# data viz
plot_annual_metric <- function(data, y_var, y_lab = y_var, plot_title = y_lab) {
  ggplot(data, aes(x = year, y = .data[[y_var]])) +
    geom_line() +
    geom_point() +
    labs(
      x = "Year",
      y = y_lab,
      title = plot_title
    ) +
    theme_bw()
}

threshold_lab <- paste0("> ", height_threshold, " cm")

predictors <- c(
  total_stems = "Total Number of Stems",
  n_plots = "Number of Plots Sampled",
  mean_stems_per_plot = "Mean Number of Stems per Plot",
  mean_height = 'Mean Height (m) of Stems Across All Plots',
  median_height = 'Median Height (m) of Stems Across All Plots',
  q90_height = '90th Quantile Height (m) of Stems Across All Plots',
  n_stems_over_threshold = paste0("Total Number of Stems ", threshold_lab),
  mean_n_stems_over_threshold_per_plot = paste0("Mean Number of Stems ", threshold_lab, " per Plot"),
  median_n_stems_over_threshold_per_plot = paste0("Median Number of Stems ", threshold_lab, " per Plot")
)

annual_plot_list <- list()

for (y_var in names(predictors)) {
  annual_plot_list[[y_var]] <- plot_annual_metric(
    data = annual_aspen_metrics,
    y_var = y_var,
    y_lab = predictors[[y_var]],
    plot_title = predictors[[y_var]]
  )
}

plot_grid(
  plotlist = annual_plot_list,
  ncol = 3
)