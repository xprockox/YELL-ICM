### ----- CLIMATE COVARIATE EXTRACTION SCRIPT FOR ICM -----
### Includes:
### -- Summer & winter precip/tmean + winter severity from PRISM data 
###    (25m spatial resolution, monthly temporal resolution)
### -- Annual NPP from MODIS product MOD17A3
###    (res)
### -- Annual brown-down / onset greenness minimum from MODIS product MCD12Q2
###    (res)
### -- Annual mean May–Aug PDSI
###    (res)
### Last updated: June 27, 2026
### xprockox@gmail.com

################################################################################
############---------------- Load packages ----------------#####################
################################################################################

library(prism)
library(terra)
library(sf)
library(lubridate)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)

################################################################################
############---------------- User settings ----------------#####################
################################################################################

prism_set_dl_dir("data/covariates/prism_data")

start_year <- 1995
end_year <- 2023

# since MODIS only goes back to 2001, it gets its own date range
modis_start_year <- 2001
modis_end_year <- end_year

elk_mcp_path <- "data/elk_gps_mcp/elk_gps_mcp_all.gpkg"

# local data directories
mod17_npp_dir <- "data/covariates/modis/MOD17A3HGF"
mcd12q2_dir <- "data/covariates/modis/MCD12Q2"
pdsi_dir <- "data/covariates/pdsi"
pdsi_path <- file.path(pdsi_dir, "noaa_pdsi_monthly.csv")

# set TRUE if you want the script to try downloading missing data automatically
download_missing_prism <- TRUE
download_missing_pdsi <- TRUE

# MODIS layer names and patterns to extract from local files
mod17_npp_pattern <- "Npp_500m"
mcd12q2_browndown_pattern <- "Dormancy_0" # alternatives: "Senescence_0", "Greenup_0"

# whether to stop if expected MODIS years are missing
stop_if_modis_incomplete <- TRUE

################################################################################
############---------------- Helper functions ----------------##################
################################################################################

# ensure directories exist
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

# ensure expected files exist
check_required_files <- function(paths) {
  tibble(
    path = paths,
    exists = file.exists(paths)
  )
}

# check to see how many of the expected files exist
check_required_count <- function(files, expected_n, label) {
  tibble(
    dataset = label,
    found = length(files),
    expected = expected_n,
    complete = length(files) >= expected_n
  )
}

# download prism
download_prism_if_missing <- function(start_year, end_year) {
  years <- start_year:end_year
  
  all_archives <- prism_archive_ls()
  
  ppt_archives <- all_archives[
    str_detect(all_archives, "prism_ppt_us_25m_") &
      str_detect(all_archives, paste0("(", paste(years, collapse = "|"), ")"))
  ]
  
  tmean_archives <- all_archives[
    str_detect(all_archives, "prism_tmean_us_25m_") &
      str_detect(all_archives, paste0("(", paste(years, collapse = "|"), ")"))
  ]
  
  expected_n <- length(years) * 12
  
  prism_status <- bind_rows(
    check_required_count(ppt_archives, expected_n, "PRISM ppt"),
    check_required_count(tmean_archives, expected_n, "PRISM tmean")
  )
  
  print(prism_status)
  
  if (length(ppt_archives) < expected_n) {
    message("Downloading missing PRISM precipitation archives...")
    for (yr in years) {
      for (mo in 1:12) {
        prism_name <- sprintf("prism_ppt_us_25m_%04d%02d", yr, mo)
        if (!prism_name %in% all_archives) {
          get_prism_monthlys(
            type = "ppt",
            years = yr,
            mon = mo,
            keepZip = FALSE
          )
        }
      }
    }
  }
  
  all_archives <- prism_archive_ls()
  
  if (length(tmean_archives) < expected_n) {
    message("Downloading missing PRISM mean temperature archives...")
    for (yr in years) {
      for (mo in 1:12) {
        prism_name <- sprintf("prism_tmean_us_25m_%04d%02d", yr, mo)
        if (!prism_name %in% all_archives) {
          get_prism_monthlys(
            type = "tmean",
            years = yr,
            mon = mo,
            keepZip = FALSE
          )
        }
      }
    }
  }
}

# download pdsi
download_pdsi_if_missing <- function(pdsi_path) {
  if (file.exists(pdsi_path)) {
    message("PDSI file already exists. Skipping download.")
    return(invisible(NULL))
  }
  
  ensure_dir(dirname(pdsi_path))
  
  # Replace this URL if you choose a different official NOAA source file
  pdsi_url <- "https://www.ncei.noaa.gov/pub/data/cirs/climdiv/climdiv-pdsidv-v1.0.0-20260604"
  
  message("Downloading NOAA PDSI file...")
  download.file(pdsi_url, destfile = pdsi_path, mode = "wb")
}

# check
get_modis_tif_inventory <- function(tif_dir, layer_pattern = NULL) {
  tif_files <- list.files(
    tif_dir,
    pattern = "\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  inv <- tibble(
    file = tif_files,
    year = as.integer(str_extract(basename(tif_files), "(19|20)\\d{2}"))
  ) %>%
    filter(!is.na(year))
  
  if (!is.null(layer_pattern)) {
    inv <- inv %>%
      filter(str_detect(basename(file), layer_pattern))
  }
  
  inv %>%
    arrange(year, file)
}

check_modis_files <- function(start_year, end_year, mod17_dir, mcd12q2_dir,
                              mod17_pattern, mcd12q2_pattern) {
  years <- start_year:end_year
  
  mod17_inv <- get_modis_tif_inventory(
    mod17_dir,
    layer_pattern = mod17_pattern
  )
  
  mcd12q2_inv <- get_modis_tif_inventory(
    mcd12q2_dir,
    layer_pattern = mcd12q2_pattern
  )
  
  tibble(
    dataset = c("MOD17A3HGF NPP", "MCD12Q2 brown-down"),
    years_expected = c(length(years), length(years)),
    years_found = c(
      sum(years %in% unique(mod17_inv$year)),
      sum(years %in% unique(mcd12q2_inv$year))
    ),
    complete = c(
      all(years %in% unique(mod17_inv$year)),
      all(years %in% unique(mcd12q2_inv$year))
    )
  )
}

################################################################################
############---------------- Data checks / downloads ---------------############
################################################################################

ensure_dir(mod17_npp_dir)
ensure_dir(mcd12q2_dir)
ensure_dir(pdsi_dir)

# check MCP exists first
mcp_check <- check_required_files(elk_mcp_path)
print(mcp_check)

if (!file.exists(elk_mcp_path)) {
  stop("Elk MCP file not found: ", elk_mcp_path)
} else {
  elk_mcp <- st_read(elk_mcp_path, quiet = TRUE) %>%
    st_make_valid()
}

# PRISM
if (download_missing_prism) {
  download_prism_if_missing(start_year, end_year)
} else {
  message("PRISM auto-download disabled. Assuming files already exist.")
}

# PDSI
if (download_missing_pdsi) {
  download_pdsi_if_missing(pdsi_path)
} else {
  message("PDSI auto-download disabled. Assuming file already exists.")
}

# MODIS local file inventory
modis_status <- check_modis_files(
  start_year = modis_start_year,
  end_year = modis_end_year,
  mod17_dir = mod17_npp_dir,
  mcd12q2_dir = mcd12q2_dir,
  mod17_pattern = mod17_npp_pattern,
  mcd12q2_pattern = mcd12q2_browndown_pattern
)

print(modis_status)

if (stop_if_modis_incomplete && any(!modis_status$complete)) {
  stop(
    paste(
      "Some MODIS years are missing from the local folders.",
      "Please make sure the GeoTIFF files for all years are present before continuing."
    )
  )
} else if (any(!modis_status$complete)) {
  warning(
    paste(
      "Some MODIS years are missing.",
      "The script will continue, but NPP / brown-down outputs may be incomplete."
    )
  )
}

# final file checks before proceeding
final_required_files <- c(
  elk_mcp_path,
  pdsi_path
)

final_file_status <- check_required_files(final_required_files)
print(final_file_status)

if (!all(final_file_status$exists)) {
  stop("One or more required non-MODIS files are still missing.")
}

################################################################################
############------------------- PRISM data -----------------------##############
################################################################################

### Find PRISM archives
all_archives <- prism_archive_ls()

ppt_archives <- all_archives[
  str_detect(all_archives, "prism_ppt_us_25m_")
]

tmean_archives <- all_archives[
  str_detect(all_archives, "prism_tmean_us_25m_")
]

### Keep only desired years
ppt_archives <- ppt_archives[
  str_detect(ppt_archives, paste0("(", paste(start_year:end_year, collapse = "|"), ")"))
]

tmean_archives <- tmean_archives[
  str_detect(tmean_archives, paste0("(", paste(start_year:end_year, collapse = "|"), ")"))
]

### Build raster stacks
ppt_stack_r <- prism_stack(ppt_archives)
tmean_stack_r <- prism_stack(tmean_archives)

ppt_stack <- rast(ppt_stack_r)
tmean_stack <- rast(tmean_stack_r)

### Align projection + crop
elk_mcp_vect <- vect(st_transform(elk_mcp, crs(ppt_stack)))

ppt_elk_mcp <- crop(ppt_stack, elk_mcp_vect) %>%
  mask(elk_mcp_vect)

tmean_elk_mcp <- crop(tmean_stack, elk_mcp_vect) %>%
  mask(elk_mcp_vect)

### Extract monthly area means
ppt_monthly <- global(ppt_elk_mcp, mean, na.rm = TRUE)[, 1]
tmean_monthly <- global(tmean_elk_mcp, mean, na.rm = TRUE)[, 1]

### Build monthly dataframe
dates <- seq(
  as.Date(sprintf("%d-01-01", start_year)),
  as.Date(sprintf("%d-12-01", end_year)),
  by = "month"
)

climate_monthly <- data.frame(
  date = dates,
  year = year(dates),
  month = month(dates),
  ppt_mm = as.numeric(ppt_monthly),
  tmean_c = as.numeric(tmean_monthly)
) %>%
  arrange(date)

print(climate_monthly)

### Build annual summaries
winter_summary <- climate_monthly %>%
  mutate(
    winter_year = ifelse(month == 12, year + 1, year)
  ) %>%
  filter(month %in% c(12, 1, 2)) %>%
  group_by(winter_year) %>%
  summarise(
    winter_ppt_mm = sum(ppt_mm, na.rm = TRUE),
    winter_tmean_c = mean(tmean_c, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(year = winter_year) %>%
  mutate(
    winter_ppt_z = as.numeric(scale(winter_ppt_mm)),
    winter_tmean_z = as.numeric(scale(winter_tmean_c)),
    winter_severity = winter_ppt_z - winter_tmean_z
  )

climate_annual <- climate_monthly %>%
  group_by(year) %>%
  summarise(
    annual_ppt_mm = sum(ppt_mm, na.rm = TRUE),
    summer_ppt_mm = sum(ppt_mm[month %in% c(6, 7, 8)], na.rm = TRUE),
    summer_tmean_c = mean(tmean_c[month %in% c(6, 7, 8)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(winter_summary, by = "year") %>%
  arrange(year)

print(climate_annual)

################################################################################
############---------------- Helper function ----------------------#############
################################################################################

extract_modis_annual_mean <- function(file_inventory, polygon_sf, years) {
  out <- lapply(years, function(yr) {
    
    th_files <- file_inventory %>%
      filter(year == yr) %>%
      pull(file)
    
    if (length(th_files) == 0) {
      return(tibble(year = yr, value = NA_real_))
    }
    
    th_rasters <- lapply(th_files, terra::rast)
    
    th_rast <- if (length(th_rasters) == 1) {
      th_rasters[[1]]
    } else {
      do.call(terra::mosaic, th_rasters)
    }
    
    poly_v <- terra::vect(sf::st_transform(polygon_sf, terra::crs(th_rast)))
    val <- terra::extract(th_rast, poly_v, fun = mean, na.rm = TRUE)
    
    tibble(
      year = yr,
      value = as.numeric(val[1, 2])
    )
  })
  
  bind_rows(out) %>%
    arrange(year)
}

################################################################################
############------------------- NPP data -------------------------##############
################################################################################

mod17_inventory <- get_modis_tif_inventory(
  mod17_npp_dir,
  layer_pattern = mod17_npp_pattern
)

npp_annual <- extract_modis_annual_mean(
  file_inventory = mod17_inventory,
  polygon_sf = elk_mcp,
  years = modis_start_year:modis_end_year
) %>%
  rename(annual_npp = value)

print(npp_annual)

################################################################################
############------------- Brown-down timing data -----------------##############
################################################################################

mcd12q2_inventory <- get_modis_tif_inventory(
  mcd12q2_dir,
  layer_pattern = mcd12q2_browndown_pattern
)

browndown_annual <- extract_modis_annual_mean(
  file_inventory = mcd12q2_inventory,
  polygon_sf = elk_mcp,
  years = modis_start_year:modis_end_year
) %>%
  rename(browndown_onset_greenness_min = value)

print(browndown_annual)

################################################################################
############----------------- Drought data -----------------------##############
################################################################################

pdsi_raw <- read.fwf(
  pdsi_path,
  widths = c(10, rep(7, 12)),
  stringsAsFactors = FALSE
)

colnames(pdsi_raw) <- c(
  "code_year",
  paste0("month_", 1:12)
)

pdsi_raw <- pdsi_raw %>%
  mutate(
    code_year = trimws(as.character(code_year)),
    climdiv_code = substr(code_year, 1, 5),
    year = as.integer(substr(code_year, 6, 10))
  )

pdsi_summer <- pdsi_raw %>%
  pivot_longer(
    cols = starts_with("month_"),
    names_to = "month",
    names_prefix = "month_",
    values_to = "pdsi"
  ) %>%
  mutate(
    month = as.integer(month),
    pdsi = as.numeric(pdsi)
  ) %>%
  filter(
    year %in% (start_year:end_year),
    month %in% 5:8
  ) %>%
  group_by(year) %>%
  summarise(
    summer_avg_pdsi = mean(pdsi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year)

################################################################################
############---------------- Combine annual outputs ---------------#############
################################################################################

environmental_covars_annual <- climate_annual %>%
  full_join(npp_annual, by = "year") %>%
  full_join(browndown_annual, by = "year") %>%
  full_join(pdsi_summer, by = "year") %>%
  arrange(year)

print(environmental_covars_annual)

################################################################################
############---------------- Export outputs ----------------------##############
################################################################################

stop("WARNING: The following lines will overwrite data. Are you sure you would like to proceed?")

# write.csv(
#   climate_monthly,
#   "data/covariates/prism_monthly_precip_tmean.csv",
#   row.names = FALSE
# )

write.csv(
  environmental_covars_annual,
  "data/covariates/environmental_covariates_annual.csv",
  row.names = FALSE
)
