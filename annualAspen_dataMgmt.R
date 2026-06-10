### Aspen Data Management for ICM
### Last updated: June 4, 2026
### xprockox@gmail.com

#################################################################################################
# packages

library(tidyverse)
library(cowplot)

#################################################################################################
# data import

dat <- read.csv('data/2024-02-25_AspenData_Clean.csv')

#################################################################################################
# data management

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
    
    n_stems_over_120 = sum(height_cm > 120, na.rm = TRUE),
    prop_over_120 = mean(height_cm > 120, na.rm = TRUE),
    mean_n_stems_over_120_per_plot = n_stems_over_120 / n_plots,
    
    mean_cag = mean(cag_cm, na.rm = TRUE),
    median_cag = median(cag_cm, na.rm = TRUE),
    
    prop_browsed = mean(browse_winter == 1, na.rm = TRUE),
    
    n_over_120_unbrowsed = sum(height_cm > 120 & browse_winter == 0, na.rm = TRUE),
    
    .groups = "drop"
  )

annual_aspen_plot_metrics <- dat %>%
  group_by(year, plot) %>%
  summarise(
    n_stems = n(),
    n_new = sum(new == 1, na.rm = TRUE),
    n_over_120 = sum(height_cm > 120, na.rm = TRUE),
    prop_over_120 = mean(height_cm > 120, na.rm = TRUE),
    median_height = median(height_cm, na.rm = TRUE),
    prop_browsed = mean(browse_winter == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(year) %>%
  summarise(
    mean_plot_stems = mean(n_stems, na.rm = TRUE),
    median_plot_stems = median(n_stems, na.rm = TRUE),
    mean_plot_over_120 = mean(n_over_120, na.rm = TRUE),
    median_plot_over_120 = median(n_over_120, na.rm = TRUE),
    prop_plots_with_over_120 = mean(n_over_120 > 0, na.rm = TRUE),
    prop_plots_with_5plus_over_120 = mean(n_over_120 >= 5, na.rm = TRUE),
    mean_plot_median_height = mean(median_height, na.rm = TRUE),
    mean_plot_browse = mean(prop_browsed, na.rm = TRUE),
    .groups = "drop"
  )


#################################################################################################
# data viz

plot_annual_metric <- function(data, y_var, y_lab, plot_title) {
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

annual_plot_info <- tibble(
  plot_name = c(
    "plot_n_stems",
    "plot_n_plots",
    "plot_mean_stems_per_plot",
    "plot_n_stems_over_120",
    "plot_mean_n_stems_over_120_per_plot"
  ),
  y_var = c(
    "total_stems",
    "n_plots",
    "mean_stems_per_plot",
    "n_stems_over_120",
    "mean_n_stems_over_120_per_plot"
  ),
  y_lab = c(
    "Total Number of Stems",
    "Number of Plots Sampled",
    "Mean Number of Stems per Plot",
    "Total Number of Stems > 120 cm",
    "Mean Number of Stems > 120 cm per Plot"
  ),
  plot_title = c(
    "Total Number of Stems",
    "Number of Plots Sampled",
    "Mean Number of Stems per Plot",
    "Total Number of Stems > 120 cm",
    "Mean Number of Stems > 120 cm per Plot"
  )
)

annual_plot_list <- annual_plot_info %>%
  mutate(
    plot = pmap(
      list(y_var, y_lab, plot_title),
      ~ plot_annual_metric(
        data = annual_aspen_metrics,
        y_var = ..1,
        y_lab = ..2,
        plot_title = ..3
      )
    )
  ) %>%
  select(plot_name, plot) %>%
  deframe()

annual_plot_list$plot_n_stems
annual_plot_list$plot_n_plots
annual_plot_list$plot_mean_stems_per_plot
annual_plot_list$plot_n_stems_over_120
annual_plot_list$plot_mean_n_stems_over_120_per_plot

plot_grid(
  plotlist = annual_plot_list,
  ncol = 3
)